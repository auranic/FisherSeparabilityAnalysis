%% Spherical Distribution Function
% This function returns an array of points that are uniformly distributed
% within the volume of a sphere of dimension dim and unit radius.
%
% The naive way to perform this work would be to first uniformly generate
% points in a cube, and then reject those outside the sphere. This becomes
% computationally costly for higher dimensional spheres, which occupy small
% fractions of their corresponding cubes. By 6D for example, only
% (pi^3/6/64) = 8% of molecules generated in a 6D-cube fall within the
% 6D-sphere.
%
% Instead, we can automatically generate points only within the hypersphere
% by working in a spherical coordinate system where this condition can be
% expressed in a single variable. It is then possible to randomly select
% each of the N-coordinates independtly, one radius and N-1 angles, with
% proper weighting achieved by inverting the desired CDF, which can be
% derived from the N-D spherical Jacobian.
%
% Points can be returned in cartesian or spherical form by setting the
% 'Coordinates' property:
% p = spheredist(1000,4,'Coordinates','Spherical')
%
% Example 1: Fill a 3D-sphere.
% p = spheredist(10000,3);
% figure;
% plot3(p(:,1),p(:,2),p(:,3),'k.')
%
% Example 2: Project a 6D-sphere into 3D. Note the increased point density
% near to the center of the observed 3D projection.
% p = spheredist(10000,6);
% figure;
% plot3(p(:,1),p(:,2),p(:,3),'k.')
%
% Example 3: Works for a simple 2D-sphere as well.
% p = spheredist(10000,2);
% figure;
% plot(p(:,1),p(:,2),'k.')
%
% Example 4: Check the volume element for a 3D-sphere
% p = spheredist(2000000,3);
% N = false(51,51,51);
% for i=1:2000000
%     sd = int8(round(p(i,:)*25.5))+26;
%     N(sd(1),sd(2),sd(3)) = true;%N(sd(1),sd(2),sd(3)) + 1;
% end
% fprintf('We find %2.2f%%, 4*pi/3/8=52.36%% actually.',100*sum(N(:))/52^3)
function xs = spheredist(N,dim,varargin)
% Get the radius firstly:
r = (rand(N,1)).^(1/dim);
% Grid for inverting CDF's:
x = linspace(0,pi,1001);
% Initialize all the angles and the cartesian coordinates. Get all the
% random values at once to prevent seeding issues.
phs = rand(N,dim-1);
xs = zeros(N,dim);
for i=1:dim-1
    % Exponent of sine term in the Jacobian. See
    % https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates
    e=dim-i-1;
    
    if e==0
        % First Angle has double domain of all others.
        phs(:,i) = phs(:,i)*2*pi;
    else
        % Integrate the Jacobian, which tells us PDF, to get the CDF.
        pint = cumsum(sin(x).^e)-sin(x).^e/2; %Indefinite Trapz of sorts
        
        % The CDF should run from 0 to 1.
        pint = pint/max(pint);
        
        if e<6
            % A uniform distribution, inverted via the CDF, gives the
            % desired distribution. See this article:
            % https://en.wikipedia.org/wiki/Inverse_transform_sampling
            phs(:,i) = interp1(pint,x,phs(:,i),'pchip');
        else
            % When e is too large, sin(x).^e starts showing several of the
            % first few points as all zero. We can no longer resolve the
            % tiny probabilities of this angle being close to 0 or pi. To
            % address this, we simply exclude the points close to 0 and pi
            % from the CDF. This only creates an error when one is sampling
            % the distribution more often than the inverse dynamic range of
            % a double array. You shouldn't notice the effect until you're
            % generating 10^16 points or so.
            mask = [true ~~diff(pint)];
            phs(:,i) = interp1(pint(mask),x(mask),phs(:,i),'pchip');
        end
    end
    
    % Multiply by the tangent of the last angle to convert it from a sin to
    % a cos. Then multiply by the cos of the new angle. This is just the
    % transformation from spherical to cartesian written in a way that can
    % be efficiently computed during the existing for-loop.
    if i>=2
        xs(:,i) = xs(:,i-1).*tan(phs(:,i-1)).*cos(phs(:,i));
    else
        xs(:,1) = r.*cos(phs(:,1));
    end
end
% The very last coordinate transformation.
xs(:,dim) = xs(:,dim-1).*tan(phs(:,dim-1));
% Handle the option for outputting shperical coordinates directly.
if nargin > 2
    if strcmp(varargin{1},'Coordinates')
        if strcmp(varagin{2},'Spherical')
            xs = [r ; phs];
        elseif ~strcmp(varargin{2},'Cartesian')
            error('Only ''Spherical'' and ''Cartesian'' coordinates supported.')
        end
    else
        error('Only ''Coordinates'' Property is supported.')
    end
end
end