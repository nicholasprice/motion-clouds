classdef motioncloud
  % This is a Matlab implementation of Laurent Perrinet's (INT - CNRS)
  % motion cloud spatio-temporal random textures.
  %
  % Motion clouds are a class of random phase textures. Here they are
  % implemented as dense mixtures of localized drifting gratings with
  % random positions.
  %
  % For a formal description see:
  %  
  %   Sanz Leon et al., Motion clouds: model-based stimulus synthesis of
  %   natural-like random textures for the study of motion perception.
  %   J. Neurophysiol. 107:3217?3226, 2012. 
  %
  % For Laurent Perrinet's Python implementation see:
  %
  %   https://github.com/NeuralEnsemble/MotionClouds.git
  %
  % Example usage:
  %
  %   m = motioncloud(256,256,120);     % 256 x 256 texels, 120 frames
  %
  %   % override default parameters
  %   m.th = pi/3;                      % mean orientation (radians)
  %   [m.Vx,m.Vy] = pol2cart(m.th,1.0); % mean horiz. and vert. speed
  %   m.sf = 32/m.Nx;                   % mean spatial frequency (32 cycles per frame)
  %   m.alpha = 1.0;                    % 1/f noise spectral density
  %   m.contrast = 0.12;                % contrast energy
  %   m.method = 'ENERGY';
  %
  %   % generate the spatio-temporal texture sequence
  %   s = m.getSequence();
  %
  %   % preview it...
  %   m.preview
  %     or 
  %   figure; colormap(gray(256));
  %   for ii = 1:m.Nt % loop over frames
  %     imagesc(s(:,:,ii)); pause(0.020);
  %     axis image
  %   end

  % Check Step 2 in getSequence - it's funky if not square. 
  
  % 2018-12-26 - Shaun L. Cloherty <s.cloherty@ieee.org>
  
  properties
    Nx = 2^8; % horiz. size (texels)
    Ny = 2^8; % vert. size (texels)

    Nt = 2^8; % frames per temporal period
    
    sf = 0.125; % mean spatial freq. (cycles per sample?)
    Bsf = 0.1; % spatial freq. bandwidth
    
    Vx = 1.0; % mean horiz. speed (spatial periods/temporal period?)
    Vy = 0.0; % mean vert. speed
    Bv = 0.5;

    th = 0.0; % mean orientation (radians)
    Bth = pi/16.0; % orientation bandwidth (radians)    
  end
  
  properties
    alpha = 0.0 % spectral envelope exponent (0.0 = white noise, 1.0 = 1/f noise)
    
    contrast = 1.0; % contrast (or contrast energy?)
    method = 'MICHELSON'; % MICHELSON or ENERGY
    
    ft0 = Inf; % spatio-temporal scaling factor (see eq. 14, p. 3221)

    logGabor = true; % use a log-Gabor kernel (instead of a standard Gabor)
    
    scaleSpace = 1; % (pixels per degree). Allows all spatial parameters to be specified in degrees. Use 1 to define sizes in pixels.
    scaleTime = 1; % (frame rate). Allows all temporal parameters to be specified in time (ms). Use 1 to define times in frames. 
  end
  
  methods
    function m = motioncloud(nx,ny,nt,varargin) % constructor
      m.Nx = nx;
      m.Ny = ny;
      m.Nt = nt;
      
      %% NEW
      % parse optional arguments...
      p = inputParser();
      p.KeepUnmatched = true;
      p.addParameter('sf', 0.125); % mean spatial freq. (cycles per pixel or cycles per degree)
      p.addParameter('Bsf', 0.1); % spatial freq. bandwidth
      
      p.addParameter('Vx', 1.0); % mean horiz. speed (spatial periods/temporal period?)
      p.addParameter('Vy', 0.0); % mean vert. speed
      p.addParameter('Bv', 0.5);
      
      p.addParameter('th', 0.0); % mean orientation (radians)
      p.addParameter('Bth',pi/16.0); % orientation bandwidth (radians)
      
%       p.addParameter('scaleSpace', 1);
      
       p.parse(varargin{:});
      args = p.Results;
      fi = fieldnames(args);
      for a=1:length(fi), m.(fi{a}) = args.(fi{a}); end
      
      
    end
    
    function env = getEnvelope(m)
      % Calculate the spatio-temporal frequency envelope.
      %
      % The envelope is the product of:
      %
      %   1. a speed envelope,
      %   2. a spatial frequency envelope,
      %   3. an orientation envelope, and
      %   4. a colour (i.e., spectral density) envelope.

      [fx,fy,ft] = m.getGrid();
      
      env = m.speedEnv(fx,fy,ft,m.Vx,m.Vy,m.Bv); % V(Vx,Vy,Bv)
      env = env .* m.freqEnv(fx,fy,ft,m.sf,m.Bsf,m.ft0,m.logGabor); % G(f0,Bf)
      env = env .* m.oriEnv(fx,fy,ft,m.th,m.Bth); % O(th,Bth)
      env = env .* m.colourEnv(fx,fy,ft,m.alpha); % C(alpha)
    end
    
    function varargout = getSequence(m,varargin)
      % Calculate spatio-temporal texture sequence.
      %
      % z = m.getSequence() returns the spatio-temporal texture sequence as
      % a 3d matrix, z (one frame per page), normalized and scaled to the
      % specified contrast.
      %
      % [~,zraw] = m.getSequence() returns the non-normalized (raw) texture
      % sequence.
      
      p = inputParser();
      p.addParameter('s',rng(), @(x) validateattributes(x,{'struct'},{'nonempty'})); % settings of the random number generator
      
      p.parse(varargin{:});
      args = p.Results;
      %
      
      nargoutchk(1,2);
      
      rng(args.s); % <-- (re-)set random number generator
        
      % 1. generate a random phase spectrum (uniform distribution)
      phase = 2*pi * rand(m.Nx,m.Ny,m.Nt);
          
      Z = exp(1j * phase);
                      
      % 2. multiply by the spatio-temporal frequency envelope
      Z = m.getEnvelope() .* Z; 

      % 3. compute the inverse fft
      Z = ifftshift(Z);
      Z(1,1,1) = 0.0; % kill the dc component
      
      z = real(ifftn(Z));

      if nargout > 1
        varargout{2} = z;
      end
      
      % 4. normalize and set contrast
      varargout{1} = m.normalize(z,m.contrast,m.method);
    end
    
    function [fx,fy,ft] = getGrid(m)
      fx = [1:m.Nx]./m.Nx - 0.5;
      fy = [1:m.Ny]./m.Ny - 0.5;
      
      ft = 0;
      if m.Nt ~= 1
        ft = [1:m.Nt]./m.Nt - 0.5;
      end
            
      [fx,fy,ft] = meshgrid(fx,fy,ft);
    end
  
    function preview(m,varargin)
        s = getSequence(m, varargin{:});

        figure
        
        colormap(gray(256));
        for ii = 1:m.Nt % loop over frames
            imagesc(s(:,:,ii));
            axis image
            pause(0.020);
        end
    end
    
    
  end % end public methods
  
  methods (Static)
    function env = speedEnv(fx,fy,ft,Vx,Vy,Bv)
      % Calculate speed envelope, V(Vx,Vy,Bv).
      %
      %   Vx - mean horizontal speed
      %   Vy - mean vertical speed
      %   Bv - horizontal and vertical speed bandwidth
      %
      % notes:
      %
      %  1. V = [Vx,Vy] = [0,1] is downward, and V = [1,0] is rightward motion.
      %  2. A speed of Vx = 1 corresponds to a mean displacement of 1/Nx
      %     per 'frame'.
      %
      %     To achieve one spatial period per temporal period, scale Vx by
      %     a factor of Nx/Nt.
      %  3. If Nt = 1, a single frame, speed V and speed bandwidth Bv are
      %     ignored.
      %  4. If the speed bandwidth Bv = 0, the fx-ft plane is used as the
      %     speed plane and you should set Vx = Vy = 0 to avoid aliasing.
      if size(ft,3) == 1 % a single frame, no motion
        env = ones(size(fx));
        return
      end
      
      if Bv == 0
        env = zeros(size(fx));
        env(ft == 0) = 1.0;
        return
      end
      
      import motionclouds.*
      
      f = motioncloud.rfreq(fx,fy,ft,[],true); % cleanDivision = true
      env = exp(-0.5*(fx*Vx + fy*Vy + ft).^2./(Bv*f).^2);
    end
    
    function env = freqEnv(fx,fy,ft,f0,Bf,ft0,logGabor)
      % Calculate radial frequency envelope, G(f0,Bf).
      %
      % f0 - mean spatial frequency
      % Bf - spatial frequency bandwidth
      import motionclouds.*
      
      if f0 == 0 || isinf(Bf)
        if logGabor
          env = motioncloud.spectEnv(fx,fy,ft,1.0,ft0); % alpha = 1.0 (1/f noise)
        else
          env = ones(size(fx));
        end
        return
      end
      
      if logGabor
        % see http://en.wikipedia.org/wiki/Log-normal_distribution
        f = motioncloud.rfreq(fx,fy,ft,[],true); % cleanDivision = true
        env = 1./f .* exp(-0.5*(log(f/f0).^2)/((log((f0+Bf)/f0)).^2));
      else
        env = exp(-.5*(motioncloud.rfreq(fx,fy,ft,ft0) - f0).^2/Bf.^2);
      end
      
    end
          
    function env = oriEnv(fx,fy,~,th,Bth)
      % Calculate the orientation envelope, O(th,Bth).
      %
      %   th - mean orientation (in radians)
      %   Bth - orientation bandwidth
      %
      % note:
      %
      % 1. Bth is the standard deviation of the Gaussian envelope which
      %    approximates the von Mises distribution for low bandwidths.
      % 2. the half-width at half height is approx. sqrt(2*Bth^2*log(2)).
    
      if isinf(Bth)
        % for large bandwidth, return a flat envelope
        env = ones(size(fx));
        return
      end
      
      if th == 0.0 && Bth == 0
        env = zeros(size(fx));
        env(fy == 0) = 1.0;
        return
      end
      
      angle = atan2(fy,fx);
      env = exp(cos(2*(angle-th))/4/Bth.^2);
    end
    
    function env = colourEnv(fx,fy,ft,alpha,ft0)
      % Calculate the colour (spectral density) envelope, C(alpha).
      %
      %   alpha = 0 white noise
      %   alpha = 1 pink (1/f) noise
      %   alpha = 2 red/brownian
      %
      % see http://en.wikipedia.org/wiki/1/f_noise
      import motionclouds.*

      if nargin < 5
        ft0 = Inf;
      end
      
      if alpha == 0.0 % white noise
        env = ones(size(fx));
      else
        f = motioncloud.rfreq(fx,fy,ft,ft0,true); % cleanDivision = true
        env = f.^(-1*alpha);
      end
    end
    
    function f = rfreq(fx,fy,ft,ft0,cleanDivision)
      % compute radial frequency... spherical coords in the
      % fourier domain
      if nargin < 5
        cleanDivision = false;
      end
      
      f2 = fx.^2 + fy.^2;

      if ~isinf(ft0)
        f2 = f2 + (ft/ft0).^2; % why?
      end
    
      if cleanDivision % avoids divide by zero...
        f2(f2 == 0) = Inf;
      end
    
      f = sqrt(f2);
    end
      
    function z = normalize(z,contrast,method)
      % Normalize to specified contrast.

      % remove mean
      z = z - mean(z(:));
      
      switch upper(method)
        case 'MICHELSON'
          z = 0.5 * z/max(abs(z(:))) * contrast + 0.5;
        case 'ENERGY'
          z = 0.5 * z/std(z(:)) * contrast + 0.5;
        otherwise
          error('Unknown method %s.',method);
      end
    end
    

  end % static methods
  

  
  
end % classdef