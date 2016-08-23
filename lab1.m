%Bi values
Bi = [0.01 0.1 1 10 100];

%index for Bi
i = 1;


%CALCULATE ROOTS
%plot the function to get a general idea of where roots are and also to confirm roots are found 
x = linspace(0,60,10000);
y = Bi(1,i) - x.*(besselj(1,x)./besselj(0,x));
plot(x,y)
axis([0 60 -4E4 8E4])
line( [0 1000], [0 0], 'Color', 'k')
grid on

%array to contain roots (zeta)
zeta = zeros(1, 60);
count = 1;

%intialise upper and lower bound of testing interval
%interval chosen must be small enough such that if the upper and lower values multiply to be greater than 0 it must mean there is no root contained
int = 0.1;
a = 0;
b = a + int;
c = (a+b)/2;

%30 roots is a assumed to be good enough
while count <= 60
    %values of function at end points and mid point
    ylow = Bi(1,i) - a*(besselj(1,a)/besselj(0,a));
    yhigh = Bi(1,i) - b*(besselj(1,b)/besselj(0,b));
    ymid = Bi(1,i) - c*(besselj(1,c)/besselj(0,c));

    %if product of end bounds of interval is negative then root must be contained
    %bisection method to find root then break
    if (ylow*yhigh) < 0        
        tol = 0.0001;
        for it = 0:100
            while abs(b-a) > tol
                c = (a+b)/2;
                ylow = Bi(1,i) - a*(besselj(1,a)/besselj(0,a));
                yhigh = Bi(1,i) - b*(besselj(1,b)/besselj(0,b));
                ymid = Bi(1,i) - c*(besselj(1,c)/besselj(0,c));
                if ymid == 0
                    zeta(1,count) = c;
                    count = count + 1;
                    break;
                elseif ( ymid*ylow < 0 )
                    b = c;                   
                else
                    a = c;
                end                
            end
            zeta(1,count) = a;
            count = count + 1;
            break;
        end
        
    %if product of end bounds of interval is 0 then one of end points must be root
    %find which end is root, move end point slightly, then break
    elseif (ylow*yhigh) == 0
        if ylow == 0
            zeta(1,count) = a;
            count = count + 1;
            a = a + 0.0001;
        else
            zeta(1,count) = b;
            count = count + 1;
            b = b + 0.0001;       
        end
    
    %if product of end bounds of interval is positive then interval does not contain root
    %do nothing
    else
    end
    
    %move interval forward
    a = b;
    b = a + int;
    c = (a+b)/2;
end

%FIXING ZETA
zeta = zeta(1:2:end);

%Fo values
Fo = [0.01, 0.1, 0.2, 0.33, 1, 2];

%mesh grid
r = -1:0.05:1;
z = linspace(0,1,41);
[R,Z] = meshgrid(r,z);


%CALCULATE TEMPERATURE
%calculate Cn values
Cn = (2./zeta).*(besselj(1,zeta)./((besselj(0,zeta).^2)+besselj(1,zeta).^2));

%calculate theta* values

%1
tempt = zeros(1,30);
theta = zeros(41);
for r = 21:41
    for n = 1:30
        tempt(1,n) = Cn(1,n)*exp(-(Fo(1,1)*(zeta(1,n)^2)))*besselj(0,(zeta(1,n)*R(1,r)));
    end
    theta(:,r) = sum(tempt);
end
for dup = 1:20
    theta(:,dup) = theta(1,42-dup);
end
    %calculate 1st order approximation
    thetaprox = zeros(41);
    for r = 21:41    
        thetaprox(:,r)= Cn(1,1)*exp(-(Fo(1,1)*(zeta(1,1)^2)))*besselj(0,(zeta(1,1)*R(1,r))); 
    end
    for dup = 1:20
        thetaprox(:,dup) = thetaprox(1,42-dup);
    end

%2
tempt2 = zeros(1,30);
theta2 = zeros(41);
for r = 21:41
    for n = 1:30
        tempt2(1,n) = Cn(1,n)*exp(-(Fo(1,2)*(zeta(1,n)^2)))*besselj(0,(zeta(1,n)*R(1,r)));
    end
    theta2(:,r) = sum(tempt2);
end
for dup = 1:20
    theta2(:,dup) = theta2(1,42-dup);
end
    %calculate 1st order approximation
    thetaprox2 = zeros(41);
    for r = 21:41    
        thetaprox2(:,r)= Cn(1,1)*exp(-(Fo(1,2)*(zeta(1,1)^2)))*besselj(0,(zeta(1,1)*R(1,r))); 
    end
    for dup = 1:20
        thetaprox2(:,dup) = thetaprox2(1,42-dup);
    end

%3
tempt3 = zeros(1,30);
theta3 = zeros(41);
for r = 21:41
    for n = 1:30
        tempt3(1,n) = Cn(1,n)*exp(-(Fo(1,3)*(zeta(1,n)^2)))*besselj(0,(zeta(1,n)*R(1,r)));
    end
    theta3(:,r) = sum(tempt3);
end
for dup = 1:20
    theta3(:,dup) = theta3(1,42-dup);
end
    %calculate 1st order approximation
    thetaprox3 = zeros(41);
    for r = 21:41    
        thetaprox3(:,r)= Cn(1,1)*exp(-(Fo(1,3)*(zeta(1,1)^2)))*besselj(0,(zeta(1,1)*R(1,r))); 
    end
    for dup = 1:20
        thetaprox3(:,dup) = thetaprox3(1,42-dup);
    end

%4
tempt4 = zeros(1,30);
theta4 = zeros(41);
for r = 21:41
    for n = 1:30
        tempt4(1,n) = Cn(1,n)*exp(-(Fo(1,4)*(zeta(1,n)^2)))*besselj(0,(zeta(1,n)*R(1,r)));
    end
    theta4(:,r) = sum(tempt4);
end
for dup = 1:20
    theta4(:,dup) = theta4(1,42-dup);
end

%5
tempt5 = zeros(1,30);
theta5 = zeros(41);
for r = 21:41
    for n = 1:30
        tempt5(1,n) = Cn(1,n)*exp(-(Fo(1,5)*(zeta(1,n)^2)))*besselj(0,(zeta(1,n)*R(1,r)));
    end
    theta5(:,r) = sum(tempt5);
end
for dup = 1:20
    theta5(:,dup) = theta5(1,42-dup);
end

%6
tempt6 = zeros(1,30);
theta6 = zeros(41);
for r = 21:41
    for n = 1:30
        tempt6(1,n) = Cn(1,n)*exp(-(Fo(1,6)*(zeta(1,n)^2)))*besselj(0,(zeta(1,n)*R(1,r)));
    end
    theta6(:,r) = sum(tempt6);
end
for dup = 1:20
    theta6(:,dup) = theta6(1,42-dup);
end



%FINAL PLOT
figure(2)

suptitle('Bi = 0.01')

subplot(3,2,1)
colormap('cool')
surf(R,Z,theta, 'FaceColor','interp', 'EdgeColor','none','FaceLighting','phong')
hold on
p = surf(R,Z,thetaprox, 'FaceColor','interp', 'EdgeColor','none','FaceLighting','phong')
set(p,'facealpha',.4,'facecolor','black')
colorbar('EastOutside')
camlight('left')
box on
axis([-1 1 0 1 0 1.1])
set(gca,'yTick',0:1:1)
set(gca,'yTickLabel',{'0','inf'})
set(gca,'xTick',-1:1:1)
set(gca,'xTickLabel',{'1','0', '1'})
title('Fo = 0.01')
xlabel('r*')
ylabel('z')
zlabel('\theta*')

subplot(3,2,2)
colormap('cool')
surf(R,Z,theta2, 'FaceColor','interp', 'EdgeColor','none','FaceLighting','phong')
hold on
p = surf(R,Z,thetaprox2, 'FaceColor','interp', 'EdgeColor','none','FaceLighting','phong')
set(p,'facealpha',.4,'facecolor','black')
colorbar('EastOutside')
camlight('left')
box on
axis([-1 1 0 1 0 1.1])
set(gca,'yTick',0:1:1)
set(gca,'yTickLabel',{'0','inf'})
set(gca,'xTick',-1:1:1)
set(gca,'xTickLabel',{'1','0', '1'})
title('Fo = 0.1')
xlabel('r*')
ylabel('z')
zlabel('\theta*')

subplot(3,2,3)
colormap('cool')
surf(R,Z,theta3, 'FaceColor','interp', 'EdgeColor','none','FaceLighting','phong')
hold on
p = surf(R,Z,thetaprox3, 'FaceColor','interp', 'EdgeColor','none','FaceLighting','phong')
set(p,'facealpha',.4,'facecolor','black')
colorbar('EastOutside')
camlight('left')
box on
axis([-1 1 0 1 0 1.1])
set(gca,'yTick',0:1:1)
set(gca,'yTickLabel',{'0','inf'})
set(gca,'xTick',-1:1:1)
set(gca,'xTickLabel',{'1','0', '1'})
title('Fo = 0.2')
xlabel('r*')
ylabel('z')
zlabel('\theta*')

subplot(3,2,4)
surf(R,Z,theta4, 'FaceColor','interp', 'EdgeColor','none','FaceLighting','phong')
colormap('cool')
colorbar('EastOutside')
camlight('left')
box on
axis([-1 1 0 1 0 1.1])
set(gca,'yTick',0:1:1)
set(gca,'yTickLabel',{'0','inf'})
set(gca,'xTick',-1:1:1)
set(gca,'xTickLabel',{'1','0', '1'})
title('Fo = 0.33')
xlabel('r*')
ylabel('z')
zlabel('\theta*')

subplot(3,2,5)
surf(R,Z,theta5, 'FaceColor','interp', 'EdgeColor','none','FaceLighting','phong')
colormap('cool')
colorbar('EastOutside')
camlight('left')
box on
axis([-1 1 0 1 0 1])
set(gca,'yTick',0:1:1)
set(gca,'yTickLabel',{'0','inf'})
set(gca,'xTick',-1:1:1)
set(gca,'xTickLabel',{'1','0', '1'})
title('Fo = 1')
xlabel('r*')
ylabel('z')
zlabel('\theta*')

subplot(3,2,6)
surf(R,Z,theta6, 'FaceColor','interp', 'EdgeColor','none','FaceLighting','phong')
colormap('cool')
colorbar('EastOutside')
camlight('left')
box on
axis([-1 1 0 1 0 1])
set(gca,'yTick',0:1:1)
set(gca,'yTickLabel',{'0','inf'})
set(gca,'xTick',-1:1:1)
set(gca,'xTickLabel',{'1','0', '1'})
title('Fo = 2')
xlabel('r*')
ylabel('z')
zlabel('\theta*')