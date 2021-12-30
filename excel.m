 %% prepare a matrix to write to a csv file
figure()

clear;

csv_matrix = zeros(1000, 9);

i = 1;
x = ["depth" "turns" "intervals/turn density" "polywidth" "capacitance (pF)" "inductance (nH)" "resistance (ohm)" "resonance (MHz)" "q_factor"];
csv_matrix = [x; csv_matrix];

i = i + 1;


% max_depth = 5;
% max_rotations = 5;
% max_intervals = 100;
% max_width = 2;
% % vary the depth of the hilbert curves from 3 - 5
% for depth = 3 : max_depth
%     %% vary the rotations of the hilbert curves from 2 - 5
%     for rotations = 2 : max_rotations
%         %% vary the intervals of the hilbert curves from 50 - 100
%         for intervals = 25 : 25 : max_intervals
%             %% vary the width between the polygon from 0.5 - 1
%             for polywidth = 0.5 : 0.25 : max_width
%                 [capacitance, inductance, resistance, resonance, q_factor] = get_measurements(false, depth, rotations, intervals, polywidth);
%                 measurements = [depth, rotations, intervals, polywidth, capacitance, inductance, resistance, resonance, q_factor];
%                 csv_matrix(i, :) = measurements;
%                 i = i + 1;
%             end
%             
%         end
%     end
% end


[capacitance, inductance, resistance, resonance, q_factor] = get_measurements(false, 3, 4, 50, 1);


max_depth = 4;
max_rotations = 9;
max_intervals = 25;
max_width = 2;

%% 3: 8 : 20 : 1.5
for depth = 3 : max_depth
    %% vary the rotations of the hilbert curves from 2 - 5
    for rotations = 7 : 0.5 : max_rotations
        %% vary the intervals of the hilbert curves from 50 - 100
        for intervals = 15 : 5 : max_intervals
            %% vary the width between the polygon from 0.5 - 1
            for polywidth = 1 : 0.2 : max_width
                [capacitance, inductance, resistance, resonance, q_factor] = get_measurements(false, depth, rotations, intervals, polywidth);
                measurements = [depth, rotations, intervals, polywidth, capacitance, inductance, resistance, resonance, q_factor];
                csv_matrix(i, :) = measurements;
                i = i + 1;
            end
            
        end
    end
end




writematrix(csv_matrix, 'hilbert_variations_1.csv')
% display(csv_matrix)
% [capacitance, inductance, resistance] = get_measurements(false, 3, 4, 50, 0.5);
% display(capacitance)
% display(inductance)
% display(resistance_spiral)



%% returns the capacitance, inductance, and resistance
%% takes in depth, number of rotations, the interval between each rotation, and the width between the polygons
function [capacitance, inductance, resistance, resonance, q_factor] = get_measurements (doPlot, depth, numRotations, interval, polywidth)
    %% Input handling:
    if nargin < 1 || isempty(doPlot)
      doPlot = false;
    end
    %% Initialization:
    origin = [0,0];
    outerD = 1000;
    thetamax = numRotations*(2*pi);
    %% Hilbert and polygon preallocate
    [hilX,hilY] = hilbert(depth);
    [xPolygon, yPolygon] = drawPolygon(hilX, hilY, polywidth);
    %currently hard-coded
    if mod(depth,2)~=0
        yPolygon(1)= yPolygon(1)+1;
        yPolygon(end)= yPolygon(end)+1;
    end
    hil_height = max(hilY);
    %preallocate for hilbert position
    width = abs(hilX(1)-hilX(end));
    DST = 2*width+2*polywidth+3; 
    %trapezoid
    %gap <1.2 micron when depth<3
    if depth>3
        shear = 2^(depth-2)-2;
        [hilX,hilY] = trapezoid(hilX,hilY,shear);
        [xPolygon,yPolygon] = trapezoid(xPolygon,yPolygon,shear);
    else
        [hilX,hilY] = trapezoid(hilX,hilY,0.3);
        [xPolygon,yPolygon] = trapezoid(xPolygon,yPolygon,0.3);
    end
    %% Calculation of (x,y) - underlying archimedean spiral. with fixed outer diameter 1000um
    innerD = outerD-2*(numRotations-1)*width-2*numRotations*(interval-width);
    innerR = outerD/2-numRotations*hil_height-numRotations*interval;
    theta = 0:0.01:thetamax; %resolution
    x = (innerR + interval.*theta/(2*pi)).*cos(theta);
    y = (innerR + interval.*theta/(2*pi)).*sin(theta);
    L = spiralL(numRotations,outerD,innerR);

    %% keep track of start and end points of hilbert polygon
    xPolygonEndPoints = [];
    yPolygonEndPoints = [];
    xPolygonStartPoints = [];
    yPolygonStartPoints = [];
    %% Calculation of equidistant (xi,yi) points on spiral.
    %length of one seg
    d = diff([hilX(:) hilY(:)]);
    hil_length = sum(sqrt(sum(d.*d,2)));

    %Start a new figure
    
    cRes = 360; 
    numSeg = ceil(2*L/DST);
    [xi,yi] = deal(NaN(numSeg,1));
    %plot(x,y,'b-','linewidth',1);
    hold on; axis equal; grid on; grid minor;
    hAx = gca; hAx.XLim = [-5 105]; hAx.YLim = [-5 105];
    xi(1) = x(end); yi(1) = y(end);

    %% calculate the final and current radius
    %final radius is equal to furthestRadius - ((numRotations-1)*interval)
    currRadius = distance(x(end), y(end));
    finalRadius = currRadius - ((numRotations-1)*interval);

    %% add hilbert curves until the abs value of the curr radius is less than or equal to the min radius
    % for ind = 2:numSeg
    ind = 2;
    line_length = 0;
    while (abs(currRadius) > abs(finalRadius))
      P = InterX([x;y], makeCircle([xi(ind-1),yi(ind-1)],DST/2,cRes));  
      absolute = abs(P(1,:)-origin(1)+1i*(P(2,:)-origin(2)));
      % get the starting point for hilbert
      [~,I] = min(absolute);
      xi(ind) = P(1,I);
      yi(ind) = P(2,I);

      % create a hilbert curve based on last intersection point
      shiftX = xi(ind-1);
      shiftY = yi(ind-1);
      [hilbertX, hilbertY] = shiftData(hilX, hilY, shiftX, shiftY);
      %find current intersection point
      curX = xi(ind);
      curY = yi(ind);
      v = [curX-hilbertX(1) curY-hilbertY(1)];

      %rotate hilbert and polygon
      ang = vec_angle(v(1), v(2));
      polyAng = vec_angle(v(1), v(2));
      [hilbertX, hilbertY] = rotxy(hilbertX, hilbertY, ang);
      [xP, yP] = rotxy(xPolygon, yPolygon, polyAng);

      % shift to proper place on spiral
      [hilbertX, hilbertY] = shiftData(hilbertX, hilbertY, shiftX, shiftY);
      [xP, yP] = shiftPolygon(xP, yP, shiftX, shiftY); 

      % add the start and end points of the curr polygon
      xPolygonEndPoints = [xPolygonEndPoints; xP(end)];
      yPolygonEndPoints = [yPolygonEndPoints; yP(end)];
      if not(ind == 2)
          xPolygonStartPoints = [xPolygonStartPoints; xP(1)];
          yPolygonStartPoints = [yPolygonStartPoints; yP(1)];
      end

      % calculate the current radius
      currRadius = avgRadius(hilbertX, hilbertY);   

      if(abs(currRadius) > abs(finalRadius))
         % draw connection lines and hilbert curve
        drawHilbert(hilbertX, hilbertY);
        plot([hilbertX(end) xi(ind)], [hilbertY(end) yi(ind)],'b','linewidth',1);
        ld = [abs(hilbertX(end)-xi(ind)) abs(hilbertY(end)-yi(ind))];
        l_length = sqrt(ld(1)^2+ld(2)^2);
        line_length = line_length+l_length;plot(xP, yP)
        plot(xP, yP,'linewidth',1);
      end
      % increment the counter
      ind = ind + 1
    end

    %plot polygon connections
    for i = 1 : length(xPolygonStartPoints)
        % plot from end to start
        plot([xPolygonEndPoints(i) xPolygonStartPoints(i)], [yPolygonEndPoints(i) yPolygonStartPoints(i)]);
    end    

    total_length = ind*hil_length+line_length;
    spiral_length = spiralL(numRotations,innerD, outerD);
    
    % calculate outer diameter & inner diameter in meters
    d_out_meters = outerD*(10^-6);
    d_in_meters = innerD*(10^-6);
    
    a = (d_in_meters + d_out_meters)/4;
    c = (d_out_meters - d_in_meters)/2;
    
    % vacuum permeability in H/m (henries per meter)
    vacuum_permeability = 4*pi*(10^-7);
    
    % n = numRotations
    
    inductance = 31.33*vacuum_permeability*(numRotations^2)*(a^2/(8*a+11*c));
    
%     dout = outerD;
%     din = dout-2*(interval-width)*numRotations-2*width*(numRotations-1);
%     u0 = 1.257*10^3; %nH/m
%     d_avg =10^-6*(dout+din)/4;
%     spiral_width =10^-6*(dout-din)/2;
% %     %https://www.translatorscafe.com/unit-converter/pt-BR/calculator/planar-coil-inductance/
%     inductance = 31.33*u0*numRotations^2*d_avg^2/(8*d_avg+11*spiral_width);
 

    
%     % resistance
%     g_res = 0.022;%microOhm-m
%     resistance = g_res*spiral_length*10^-6/((interval-width)*1*10^-12);
%     resistance = resistance* 1*10^-6;


    resistivity = 0.022;%microOhm/m
    length_meters = spiral_length*10^-6; % convert micrometers to meters
    width_meter = (interval + width)*10^-6; % get the width of a spiral iteration and convert micrometers to meters
    thickness = 1*10^-6; % a thickness of 1 micrometer
    resistance = resistivity*length_meters/(width_meter*thickness);
    
   

%     resistivity = 0.022;
%     resistance_one = resistivity * (spiral_length) * (1/((dout+width)*10^-12));
%     
%     resistance_half = resistivity * (spiral_length*10^-6) * (1/((dout+width)*1/2*10^-12));
    
    %capacitance
    capacitance = 8.854*10^-12*polywidth*total_length*10^-6*1*10^-6/(1*10^-6); % (F)

    
    % resonance
    resonance = 1/(2*pi*sqrt(inductance*capacitance));
    
    % q-factor
%     q_factor = 1/R*sqrt(inductance/capacitance);
    q_factor = 1/resistance*sqrt(inductance/capacitance);
    
    
    % convert to ohms
    resistance = resistance * 10^-6;
    
    % convert to pF
    capacitance = capacitance * 10^12;
    
    % convert to nH
    inductance = inductance * 10^9;
    
    
    
    
%     q_factor_one = 1/(resistance_one) * sqrt(inductance/capacitance);
%     q_factor_half = 1/(resistance_half) * sqrt(inductance/capacitance);

    % disp(endx)
    % for i = 1:idx-2
    %     plot([endx(1) startx(idx+1)],[endy(1) starty(idx+1)],'r');
    %     pause(0.01);
    % end

    % delete previous Hilbert curve - went over
%     delete(hAx.Children(1))
    xi = xi(~isnan(xi)); yi = yi(~isnan(yi));



%% Nested function(s):
function [XY] = makeCircle(cnt, R, nPts)
  P = (cnt(1)+1i*cnt(2))+R*exp(linspace(0,1,nPts)*pi*2i);
  if doPlot, plot(P,'Color',lines(1)); end
  XY = [real(P); imag(P)];
end

function L = spiralL(turns, innerD, outerD)%spiral length
    L = pi*turns*(innerD+outerD)/2;
end

%% Hilbert curve:
function [XY] = drawHilbert(x,y)
  line(x, y,'linewidth',1);
  XY = [x y];
end


function [x,y] = hilbert(depth)
    A = zeros(0,2);
    B = zeros(0,2);
    C = zeros(0,2);
    D = zeros(0,2);

    north = [ 0  2];
    east  = [ 2  0];
    south = [ 0 -2];
    west  = [-2  0];

    order = depth;
    for n = 1:order
      AA = [B ; north ; A ; east  ; A ; south ; C];
      BB = [A ; east  ; B ; north ; B ; west  ; D];
      CC = [D ; west  ; C ; south ; C ; east  ; A];
      DD = [C ; south ; D ; west  ; D ; north ; B];
      A = AA;
      B = BB;
      C = CC;
      D = DD;
    end
    
    A = [0 0; cumsum(A)];
    
    xHilbert= A(:,1);
    yHilbert= A(:,2);
    x = xHilbert;
    y = yHilbert;
  
end

end

%% Local function(s):
%% funtion to find the length of the spiral
function L = spiralL(turns, innerD, outerD)%spiral length
    L = pi*turns*(innerD+outerD)/2;
end
%% function to draw a polgyon
function [xPlots, yPlots] = drawPolygon(x, y, distance)
% get the start
c1 = [x(1) y(1)];
c2 = [x(2) y(2)];
dir = getDirection(c1(1), c2(1), c1(2), c2(2));
if(strcmp(dir, "up"))
   xStart = x(1) - distance;
   xPlots = [xStart];
   yPlots = [0];
else
    yStart = y(1) + distance;
    xPlots = [0];
    yPlots = [yStart];
end

counter = 2;
disp(size(x))
lengthOfArray = size(x);
while(counter < lengthOfArray(1))
    c1 = [x(counter-1) y(counter-1)];
    c2 = [x(counter) y(counter)];
    c3 = [0 0];
    c4 = [0 0];
    c5 = [0 0];
    
    if(counter < lengthOfArray(1) - 2)
        c3 = [x(counter + 1) y(counter + 1)];
        c4 = [x(counter + 2) y(counter + 2)];
        c5 = [x(counter + 3) y(counter + 3)];
    elseif(counter < lengthOfArray(1) - 1)
        c3 = [x(counter + 1) y(counter + 1)];
        c4 = [x(counter + 2) y(counter + 2)];
    elseif(counter < lengthOfArray(1))
        c3 = [x(counter + 1) y(counter + 1)];
    end
    
    % x-components
    c1X = c1(1);
    c2X = c2(1);
    c3X = c3(1);
    c4X = c4(1);
    c5X = c5(1);
    
    % y-components
    c1Y = c1(2);
    c2Y = c2(2);
    c3Y = c3(2);
    c4Y = c4(2);
    c5Y = c5(2);
    
    % first, find the direction
    dir = getDirection(c1X, c2X, c1Y, c2Y);
    
    changeLength = determineLength(dir, c2X, c3X, c4X, c5X, c2Y, c3Y, c4Y, c5Y, distance);
    
    % get the next coordinate and add to the list of plots
    [nextX, nextY] = plotNextCoord(dir, xPlots, yPlots, c3, changeLength);
    xPlots = [xPlots nextX];
    yPlots = [yPlots nextY];
    hold on
    % if movement = up/down, add 2 + changeLength to prev X
    counter = counter + 1;
    
end

c1 = [x(end-1) y(end-1)];
c2 = [x(end) y(end)];
dir = getDirection(c1(1), c2(1), c1(2), c2(2));
if(strcmp(dir, "down"))
   xEnd = x(end) + distance;
   yEnd = y(end);
else
    xEnd = x(end);
    yEnd = y(end) + distance;
end

xPlots = [xPlots xEnd];
yPlots = [yPlots yEnd];
%  plot(xPlots, yPlots)


end

function dir = getDirection(c1X, c2X, c1Y, c2Y)
    % up
    if(c1Y < c2Y)
        dir = "up";        
    end
    % down
    if(c1Y > c2Y)
       dir = "down"; 
    end
    % left
    if (c1X > c2X)
        dir = "left";
    end
    % right
    if(c1X < c2X)
        dir = "right";
    end
    
end

function changeLength = determineLength(direction, c2X, c3X, c4X, c5X, c2Y, c3Y, c4Y, c5Y, distance)
    % assuming a left/right direction:
    if(strcmp(direction, "left") || strcmp(direction, "right"))
       % hit a wall if c3Y > c2Y - use a shorter length
       if(c3Y > c2Y)
           changeLength = -distance;
       % pass over a wall if c3Y < c2Y - use a longer length
        elseif(c3Y < c2Y)
           changeLength = distance;
       % look at the next, next coordinates if this is a longer side
       % hit a wall
        elseif(c4Y > c2Y)
           changeLength = -distance;
       % pass over a wall
        elseif(c4Y < c2Y)
           changeLength = distance;
       % special case with the longest side
       % hit a wall
       elseif(c5Y > c2Y)
           changeLength = -distance;
       % pass over a wall
       elseif(c5Y < c2Y)
           changeLength = distance;
       else
            changeLength = 0;
       end
    end
    % assuming a down direction
    if(strcmp(direction, "up") || strcmp(direction, "down"))
        % hit a wall if c2X > c3X - use a shorter length
        if(c2X > c3X)
            changeLength = -distance;
        % pass over a wall if c2X < c3X - use a longer length
        elseif(c2X < c3X)
            changeLength = distance;
        % look at the next, next coordinate if this is a longer side
        % hit a wall if c2X > c4X - use a shorter length
        elseif(c2X > c4X)
            changeLength = -distance;
        % pass over a wall if c2X < c4X - use a longer length
        elseif(c2X < c4X)
            changeLength = distance;
        % special case with the longest side
        % hit a wall if c2X > c5X - use a shorter length
        elseif(c2X > c5X)
            changeLength = -distance;
        % pass over a wall if c2X < c5X - use a longer length
        elseif(c2X < c5X)
            changeLength = distance;
        else
            changeLength = 0;
        end
    end
end


function [xCoord,yCoord] = plotNextCoord(dir, xPlots, yPlots, c3, changeLength)
    % get the last coordinates
    xLast = xPlots(end);
    yLast = yPlots(end);    
    
    xLength = c3(1) - xLast;
    yLength = c3(2) - yLast;
    
    plotLength = [0, 0];
        
    if(strcmp(dir, "left") || strcmp(dir, "right"))
        plotLength = [xLength+changeLength, 0];
    end
    if(strcmp(dir, "up") || strcmp(dir, "down"))
        plotLength = [0, yLength+changeLength];
    end
    
    
    
    % add the changes
    xCoord = xLast + plotLength(1);
    yCoord = yLast + plotLength(2);
    
       
    
end
%% function for hilbert shift 


function [shiftedX, shiftedY] = shiftData(dataX, dataY, shiftX, shiftY)
  xIncrementation = (shiftX-dataX(1));
  yIncrementation = (shiftY-dataY(1));
  shiftedX = dataX+xIncrementation;
  shiftedY = dataY+yIncrementation;
end
function [shiftedX, shiftedY] = shiftPolygon(xPolygon, yPolygon, shiftX, shiftY)
    shiftedX = xPolygon + shiftX;
    shiftedY = yPolygon + shiftY;
end
%% Fuctions to calculate the radius
function dist = distance(x, y)
    dist = sqrt(x^2 + y^2);
end

function avg = avgRadius(hilbertX, hilbertY)
    xStart = hilbertX(1);
    yStart = hilbertY(1);
    
    xEnd = hilbertX(end);
    yEnd = hilbertY(end);
    
    startRadius = distance(xStart, yStart);
    endRadius = distance(xEnd, yEnd);
    
    avg = (startRadius + endRadius)/2 ;
end

%% Functions to rotate the Hilbert curves

function [xr, yr] = rotxy(x, y, ang)
    xr = x*cos(ang)-y*sin(ang);
    yr = x*sin(ang)+y*cos(ang);
end

function alpha = vec_angle(x, y)
   if(x == 0)
      if(y >= 0)
         alpha= pi/2; 
         
      else
         alpha = -pi/2;
      end
   else if(x < 0)%case when vector in 1,3 quadrant    
    alpha = (pi)+(atan(y/x));
   else
    alpha = atan(y/x);   
       end
   end
   %alpha = alpha - pi/2;%rotate the starting-ending edge of segment
end
%% Function for making the curve and polygon trapezoid
function [trapX,trapY] = trapezoid(hx,hy,d)
height = max(hy);
midx = max(hx)/2;
hxy(:,1) = hx; hxy(:,2) = hy;
for i = 1:size(hxy,1)
    %base displacement determined by y position
    d_base = d*hxy(i,2)/height;
    %adjust based on x position
    dx = d_base*(abs(hxy(i,1)-midx)/midx);
    if hxy(i,1)>midx
        hxy(i,1) = hxy(i,1)+dx;
    end
    if  hxy(i,1)<midx
        hxy(i,1) = hxy(i,1)-dx;
    end
end
trapX = hxy(:,1); trapY = hxy(:,2);
end
%% Find the intersection between curr hilbert and spiral curve
function P = InterX(L1,varargin)
    % DOCUMENTATION REMOVED. For a full version go to:
    % https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections

    narginchk(1,2);
    if nargin == 1
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end

    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);

    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);

    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i), P = zeros(2,0); return; end

    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0

    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';

    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end