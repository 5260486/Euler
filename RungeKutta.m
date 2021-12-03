%% 龙格库塔

format long;
a=0;
b=3000;     %仿真范围
h=0.0001; %步长
x0=2; %方程初值
y0=0;
h1=0.001;
h2=0.0001;
tic
% result=Exchanged(a,b,h,x0,y0);
result=Multirate(a,b,h1,h2,x0,y0);
toc
t=result.H;
x=result.X;
plot(t,x,color);
hold on;%在同一张图绘制曲线

[t,x2]=ode45('func',[0 b],[2;0]);
plot(t,x2(:,1),'k');%精确解
legend('dt=0.0001','dt=0.0002','True');

%% 需求解常微分方程
function dx=func1(x,y)     
    dx=y;
end

function dy=func2(x,y)     
    dy=-1000*(x^2-1)*y-x;
end


%% 龙哥库塔，二阶/改进欧拉
function result=Runge2(a,b,h,x0,y0)
    n=(b-a)/h;
    X=x0;
    Y=y0;
    H=0;
        for i=1:n
            k11=h*func1(x0,y0);
            k12=h*func2(x0,y0);
            k21=h*func1(x0+h,y0+k12);
            k22=h*func2(x0+h,y0+k12);
            x1=x0+0.5*(k11+k21);
            y1=y0+0.5*(k12+k22);

            x0=x1;
            y0=y1;
            X(end+1)=x0;
            Y(end+1)=y0;
            H(end+1)=h*i;
        end 
    result.X=X;
    result.Y=Y;
    result.H=H;
end

function result=Multirate(a,b,h1,h2,x0,y0)
    n1=(b-a)/h1;%h1=0.001慢变 h2=0.0001快变;
    n2=(b-a)/h2;
    d=n2/n1;%插值位置
    X=x0;
    Y=y0;
    H=0;
        for i=1:n1           %x变化

            for j=1:d                %y变化
                k12=h2*func2(x0,y0);                          
                k22=h2*func2(x0,y0+k12);               
                y1=y0+0.5*(k12+k22);
                y0=y1;
                Y(end+1)=y0;
            end
            k11=h1*func1(x0,y0);            
            k21=h1*func1(x0+h1,y0);   
            x1=x0+0.5*(k11+k21);
            x0=x1;          
            X(end+1)=x0;           
            H(end+1)=h1*i;
        end 
    result.X=X;
    result.Y=Y;
    result.H=H;
end

%% Runge-Kutta法，四阶
function result=Runge4(a,b,h,x0,y0)
    n=(b-a)/h;
    X=x0;
    Y=y0;
    H=0;
        for i=1:n
            k11=h*func1(x0,y0);
            k12=h*func2(x0,y0);
            k21=h*func1(x0+h/2,y0+h*k12/2);
            k22=h*func2(x0+h/2,y0+h*k12/2);
            k31=h*func1(x0+h/2,y0+h*k22/2);
            k32=h*func2(x0+h/2,y0+h*k22/2);
            k41=h*func1(x0+h,y0+h*k32);
            k42=h*func2(x0+h,y0+h*k32);
 
            x1=x0+(k11+2*k21+2*k31+k41)/6;
            y1=y0+(k12+2*k22+2*k32+k42)/6;

            x0=x1;
            y0=y1;
            X(end+1)=x0;
            Y(end+1)=y0;
            H(end+1)=h*i;
        end 
    result.X=X;
    result.Y=Y;
    result.H=H;
end
