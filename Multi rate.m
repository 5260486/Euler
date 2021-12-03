format long;
a=0;
b=3000;     %仿真范围
x0=2; %方程初值
y0=0;

h1=0.001;%慢变步长
h2=0.0001;%快变步长
result=IBMEuler(a,b,x0,y0,h1,h2);
x=result.X;
y=result.Y;
t1=result.H1;
t2=result.H2;
plot(t1,x,'color','red');
hold on;

[t,x2]=ode45('func',[0 b],[2;0]);
plot(t,x2(:,1),'k');%精确解

%% 需求解常微分方程
function dx=func1(x,y)
    dx=y;
end

function dy=func2(x,y)
    dy=-1000*(x^2-1)*y-x;
end

%% 以个体为本的后向欧拉法
function result=IBMEuler(a,b,x0,y0,h1,h2)
    n1=(b-a)/h1;%h1=0.001慢变 h2=0.0001快变;
    n2=(b-a)/h2;
    d=n2/n1;%插值位置
    X=x0;
    Y=y0;
    H1=0;
    H2=0;
        for i=1:n1  %快变量变化d次，慢变量再变化           
            for j=1:d
                fy=calculate1(x0,y0,h2);
                y0=fy;            
                Y(end+1)=y0;
            end
            fx=calculate2(x0,y0,h1);
            x0=fx; 
            X(end+1)=x0;
            H1(end+1)=h1*i;
        end 
        for k=1:n2
            H2(end+1)=h2*k;
        end
    result.X=X;
    result.Y=Y;
    result.H1=H1;
    result.H2=H2;
end

function result = calculate1(x0,y0,h2)
    acc = 6;
    now = 0;      
    iter=0;
    y1=y0+h2*func2(x0,y0);
    while now <=acc
        y1=y0+h2*func2(x0,y1);
        y0=y1;
        now =abs(log10(abs(y1-y0)));   
        iter=iter+1;
        if iter>5
            break;
        end
    end         
    result = y0;
end

function result = calculate2(x0,y0,h1)
    acc = 6;
    now = 0;     
    iter=0;
    x1=x0+h1*func1(x0,y0);
    while now <=acc
        x1=x0+h1*func1(x1,y0);
        x0=x1;
        now =abs(log10(abs(x1-x0)));   
        iter=iter+1;
        if iter>5
            break;
        end
    end         
    result = x0;
end
