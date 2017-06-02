classdef BackProjection2D
    %BACKPROJECTION2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        angles  % number of angles sampled
        si      % number of sample points
        projections % projections
        BPI     % image reconstructed by back projection
        range   % range used by back projection
        
        collection % a collection of intermidiate results for a movie
    end
    
    methods
        function obj = BackProjection2D(projections, range)
            [obj.angles, obj.si] = size(projections);
            obj.projections = projections;
            obj.range = range;
            obj.si = length( obj.range );

            obj = obj.BPrecon(obj.range);
            
            obj.makeMovie();
            obj.showImage();

        end
        
        function obj = BPrecon(obj, range)
            scope = 180.0 - 180.0 / obj.angles;
            obj.range = range;
            
            PR = obj.projections(:, obj.range)';
            
            % "a":	This parameter varies the filter magnitude response.
            %	When "a" is very small (a<<1), the response approximates |w|
            %	As "a" is increased, the filter response starts to
            %	roll off at high frequencies.
            
            % a was 20
            a = 20;

            w = pi * ( -1 : (2/obj.si) : (1- 2/obj.si) )';
            
            rd = (a/2) * w;
            rn2 = sin( rd );
            rn1 = abs( (2 / a) * rn2 );
            r = rn1 .* (rn2 ./ rd).^2;
            
            f = fftshift(r);
            for k = 1:obj.angles
                IMG = fft(PR(:,k));
                fimg = IMG .* f;
                filtPR(:,k) = ifft(fimg);
            end
            filtPR = real(filtPR);
            
            th = (scope*pi/180/(obj.angles-1))*( 0:(obj.angles-1) ); % in radians
            
            BPI = zeros(obj.si, obj.si);
            
            midindex = (obj.si + 1) / 2;
            
            x = 1:obj.si;
            y = 1:obj.si;
            [X, Y] = meshgrid(x,y);
            xpr = X - midindex;
            ypr = Y - midindex;
            
            clear X Y x y
            
            collection = zeros( obj.si, obj.si, obj.angles);
            
            for k = 1:obj.angles
                disp(['Backprojecting angle ', num2str(360/(2*pi)*th(k))]);
                BPIa = zeros(obj.si,obj.si);
                
                filtIndex = round( midindex + xpr * sin(th(k)) - ypr * cos(th(k)) );
                spota = find( (filtIndex > 0) & (filtIndex <= obj.si) );
                newfiltIndex = filtIndex( spota );
                
                BPIa(spota) = PR(newfiltIndex(:),k);
                collection(:,:,k) = BPIa;
            end
            
            obj.collection = collection;
            obj.BPI = sum(collection, 3);

        end
        
        function frames = makeMovie(obj)
            frames(obj.angles) = struct('cdata',[],'colormap',[]);
            
            BPI = zeros(obj.si, obj.si);
            
            figure
            for k = 1:obj.angles
                BPI = BPI + obj.collection(:,:,k);
                imagesc(BPI)
                axis image
                axis off
                frames(k) = getframe;
            end
        end
        
        function showImage(obj)                         
            figure
            imagesc(abs(obj.BPI))
            axis image
            axis off
        end
        
    end
    
end

