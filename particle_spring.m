clear all;
close all;
axis equal;
axis tight manual;
axis([0 20 0 10]);

%% -- Static Variables --

% Variables for making the structure
XDIM   = 1;
YDIM   = 2;
dims   = [7 3];
amount = [8 4];

% Variables that control the simulation
mass = 1;
ks   = 100;   % Spring constant.
kd   = 5;     % Spring damping.
g    = -1;    % Gravitational acceleration.
dt   = 0.002; % Timestep for the simulation.
loopEnd = 3/dt;
% Rest length of springs - distance between connected particles if left 0
overrideSpringLength = 0; 

% Initial values
% Initial velocity 
iv      = zeros(prod(amount), 2);
%iv(1,:) = [0,-5];
%iv(2,:) = [0,5];
iv(:,1) = 7;
% Hight of the lower part of the object. 
h       = 0.5;

% Variables of the balls
cr           = 3/4; % Radius of the balls
randomFactor = 0;   % Multiplied by randn and added to the radius.
circleAmount = 32;

% For animating
filename = "pulse";
animate  = false;
graphs   = [1 1 0]; 
pR       = 0.03;    % Raduis for the balls representing the particle masses;
draw     = true;
ddt      = 10;      % Draws every ddt time steps.




%% -- Make circles --
circles = zeros(1,circleAmount);
cX = zeros(circleAmount, 2); % Position. 
cR = rand([1, circleAmount])*randomFactor + cr;
for i = 1:circleAmount
    cX(i,:) = [0.8*cr/2+(cr/2*(i-1)*2)+0.1*cr/2*randn() cR(i)];
    circles(i) = rectangle('Position',...
            [cX(i,XDIM) - cR(i),...
             cX(i,YDIM) - cR(i),...
             cR(i)*2, cR(i)*2], ...
            'Curvature', [1,1], 'EdgeColor','b', 'FaceColor','w');
end

%% -- Initialize particles and springs --
division = amount-1;
division(division < 1) = 1;
lDelta    = dims./division;
pX = zeros(prod(amount), 2); % Particle position vector
pV = zeros(prod(amount), 2); % Particle velocity vector
pD = zeros(prod(amount), 2); % Particle last distance to nearest floor circle
particles = zeros(prod(amount));
springs   = struct('fromTo',{}, 'line'    ,{}, 'len', {});

% Bulding particles
for x = 1:1:amount(1)
    for y = 1:1:amount(2)
        thisPoint = (y-1)*amount(1) + x;
        pX(thisPoint,:) = [((x-1) *lDelta(1)) ((y+h)*lDelta(2))];
        pV(thisPoint,:) = iv(thisPoint,:);
    end
end

% Modefiers for spring locations
width = amount(1);
sumParts = prod(amount);
mods = [-width-1, width-1, -width, width, width+1, +1, 1-width];
spring = 1;

% Particles on odd rows connect to their right neighbor
% O O O
% O P -
% O O O
%
% Particles on even rows connect to all neighbors except teir left.
% \ | /
% O P -
% / | \

% Building springs
for i = 1:sumParts
    % Checks if the row is even or odd
    if (mod(((i-1)-mod((i-1),width))/width, 2) ~= 0)
        %Filter edgecases
        usedMods = mods;
        if (mod(i,width) == 1) % Filters out the 2 leftmost particles
            usedMods = mods(3:end);
        elseif (mod(i,width) == 0)
            usedMods = mods(1:end-3); %Filters out the 3 rightmost particles
        end
        
        % Creates all the springs with their from and to particle indecies.
        for il = 1:size(usedMods,2)
            if(i + usedMods(il) > 0 && i + usedMods(il) <= sumParts)
                springs(spring).fromTo = [i, (i + usedMods(il))];
                springs(spring).force = 0;
                if overrideSpringLength == 0
                    % Sets the length to the distance bewteen the particles
                    l = norm((pX(i,:)-pX(i+usedMods(il),:)));
                    springs(spring).len = l;
                else
                    springs(spring).len = overrideSpringLength;
                end
                spring = spring + 1;
            end
        end

    else
        % Checks so that the particle is not on the far right
        if ~(mod((i-1),width) == width-1)
            springs(spring).fromTo = [i, (i + 1)];
            if overrideSpringLength == 0
                l = norm((pX(i,:)-pX(i+1,:)));
                springs(spring).len = l;
            else
                springs(spring).len = overrideSpringLength;
            end
            spring = spring + 1;
        end
    end
end

% Make graphics objects
parts   = size(particles,2);
springc = size(springs,2);

for i = 1:springc
    fromTo = springs(i).fromTo;
    springs(i).line  = line([pX(fromTo(1),XDIM), pX(fromTo(2),XDIM)],...
        [pX(fromTo(1),YDIM),pX(fromTo(2),YDIM)]);
end

for i = 1:parts
    particles(i) = rectangle('Position',...
        [pX(i,XDIM)-pR, pX(i,YDIM)-pR, 2*pR, 2*pR], ...
        'Curvature', [1,1], 'EdgeColor','r', 'FaceColor','r');
end

%% -- Initial half step back --

force = zeros (size(particles,2), 2);

% Calculate spring forces
for i = 1 : size(springs,2)
    from = springs(i).fromTo(1);
    to   = springs(i).fromTo(2);
    
    r = pX(from,:) - pX(to,:);
    v = pV(from,:) - pV(to,:);
    d = norm(r);

    f = -(ks*(d-springs(i).len) + kd*(sum(r.*v)/d))*r./d;

    force(from,:) = force(from,:) + f;
    force(to,:)   = force(to,:)   - f;
end

% Add gravity
force(:,2) = force(:,2) + g;

pV   = (pV - force/mass*dt*1/2);
for i = 1 : size(particles,2)
    if pD(i) < 2*cr
        for ci = 1 : size(cX,1)
            r = pX(i,:) - cX(ci,:);
            d = norm(r);
            pD(i) = d;

            if(d <= cR(ci))
                pX(i,:) = pX(i,:) - 1.1*(cR(ci) - d)*r/d;
                pV(i,:) = pV(i,:) - 2*(sum(pV(i,:).*r)*r/d);
            end
        end
    else
        pD(i) = 0;
    end
end

%% -- Main Loop --

% All the energy vectors
eS = ceil(loopEnd/ddt);
eK = zeros(eS,1);
eP = zeros(eS,1);
mC = zeros(eS,1);
aM = zeros(eS,1);
mI = zeros(eS,1);
iL = zeros(eS,1);
omega = zeros(eS,1);
eS = zeros(eS,1);
h = figure(1);  %For saving the frames for animation

for t = 1:1:loopEnd
    
    %% -- Update posistions and velocities --
    force = zeros (size(particles,2), 2);
    
    for i = 1 : size(springs,2)
        from = springs(i).fromTo(1);
        to   = springs(i).fromTo(2);
        
        r = pX(from,:) - pX(to,:);
        v = pV(from,:) - pV(to,:);
        d = norm(r);
        
        if(mod(t,ddt) == 0)
            eS(t/ddt) = eS(t/ddt) + 0.5 * ks * (springs(i).len - d)^2; 
        end
        
        f = -(ks*(d-springs(i).len) + kd*(sum(r.*v)/d))*r./d;

        force(from,:) = force(from,:) + f;
        force(to,:)   = force(to,:)   - f;
    end

    force(:,2) = force(:,2) + g;

    pV = (pV + force/mass*dt*1/2);
    % Saves energies between halfsteps so that I don't need to fix it when
    % saving / displaying.
    if (mod(t,ddt) == 0)
        eK(t/ddt) = sum(vecnorm(pV,2,2).^2*mass/2);
        % Saves the angular momentum stuff for the second asignment.
        if dims(2) == 0
            aM(t/ddt) = sum(mass*(pX(:,XDIM).*pV(:,YDIM) - pX(:,YDIM).*pV(:,XDIM)));
            mI(t/ddt) = sum(vecnorm([0.9 0; 0.9 0] - pX,2,2)*mass);
            omega(t/ddt) = aM(t/ddt)/mI(t/ddt);
            iL(t/ddt) = sum(vecnorm([0.9 0; 0.9 0] - pX,2,2));
        end
        mC(t/ddt) = sum(pV(:,1)*mass,1)/(mass*prod(amount));
    end
    pV = (pV + force/mass*dt*1/2);
    pX = (pX + pV*dt);
    
    %% -- Check for colisions --  
    for i = 1 : size(particles,2)
        % Only if the particle was closer than 2 radia to the nearest ball
        % last update
        if pD(i) < 2*cr
            for ci = 1 : size(cX,1)
                r = pX(i,:) - cX(ci,:);
                d = norm(r);
                
                if (pD(i)>d) 
                    pD(i) = d; 
                end
                
                if(d < cR(ci))
                    v = norm(pV(i,:));
                    nv = pV(i,:)/v;
                    % Do while that moves the particle backwards to just 
                    % outside the ball in the direction it came from.
                    % Quickfix that could take far to long depending on
                    % amount of collisions and how far into the objects the
                    % particle moves as well as how small the objects are.
                    while 1 
                        pX(i,:) = pX(i,:) - 0.01*nv;
                        r = pX(i,:) - cX(ci,:);
                        d = norm(r);
                        if ~(d <= cR(ci))
                            break
                        end
                    end
               
                    pV(i,:) = pV(i,:) - 2*(sum(pV(i,:).*r/d)*r/d);

                end
            end
        else
            % Resets the distance counter every other update.
            pD(i) = 0;
        end
    end
    
    %% -- Draw the changes --
    if(mod(t,ddt) == 0)
        
        if(draw)
            parts   = size(particles,2);
            springc = size(springs,2);

            for i = 1:parts
                p = particles(i);
                set(p, 'Position',...
                    [pX(i,XDIM)-pR, pX(i,YDIM)-pR, 2*pR, 2*pR], ...
                    'Curvature', [1,1], 'EdgeColor','r', 'FaceColor','r');
            end

            for i = 1:springc
                fromTo = springs(i).fromTo;
                set(springs(i).line,...
                    'XData', [pX(fromTo(1),XDIM), pX(fromTo(2),XDIM)],...
                    'YData', [pX(fromTo(1),YDIM), pX(fromTo(2),YDIM)]);
            end
            
            
            pause(0.016)
        end
        if animate
            drawnow
            frame(t/ddt) = getframe(h);
        end
    end
end

%% -- Graphing --
xax = 1:ddt:loopEnd;
if graphs(1)
    figure
    plot(xax', eK);
    hold on;
    plot(xax', eP);
    plot(xax', eS);
    eT = eP + eK + eS;
    plot(xax', eT);
    plot(xax', iL);
    hold off;
    legend('eK','eP','eS','eT', 'length');
end
if graphs(2)
    figure
    plot(xax', mC);
    mcv = mC(1) - mC(end)
    loopEnd*dt
    mcv = mcv/(loopEnd*dt)
end
if graphs(3)
    figure
    plot(xax', aM);
    hold on;
    plot(xax', mI);
    plot(xax', omega);
    plot(xax', iL);
    legend('aM','mI','omega', 'length');
end


if animate
    gifWriter = VideoWriter(filename);
    gifWriter.FrameRate = 30;
    open(gifWriter);
    writeVideo(gifWriter, frame);
    close(gifWriter);
end