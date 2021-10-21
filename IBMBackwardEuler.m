format long;
a=0;
b=3000;     %仿真范围
h=0.0001; %步长
x0=2; %方程初值
y0=0;

% for i=1:2
% h=h*i;%变步长
% result1=Euler(a,b,h,x0,y0);
% x=result1(1,:);
% % y=result1(2,:);
% t=result1(3,:);
%     if i==1
%         color='blue';
%     else
%         color='red';
%     end
% plot(t,x,color);
% hold on;%在同一张图绘制曲线
% end

h1=0.001;%慢变步长
h2=0.0001;%快变步长
result=IBMEuler(a,b,x0,y0,h1,h2);
x=result.X;
y=result.Y;
t1=result.H1;
t2=result.H2;
plot(t1,x,'color','green');
hold on;

[t,x2]=ode45('func',[0 b],[2;0]);
plot(t,x2(:,1),'k');%精确解

%% 需求解常微分方程
function fxy=func(x,y)  
    fxy=[y(1);
        -1000*(x(1)^2-1)*y(1)-x(1)];
end

%% 后向欧拉法
function result=Euler(a,b,h,x0,y0)
    n=(b-a)/h;
    result(1,1)=x0;
    result(2,1)=y0;
    result(3,1)=h;
    f0=func(x0,y0);
        for i=1:n
            x1=x0+h;
            f1=[x0+h*f0(1),y0+h*f0(2)];
            fy=calculate(x0,y0, x1, f1, h);
            
            f0=fy;
            x0=fy(1);
            y0=fy(2);
            
            result(1,i+1)=x0;
            result(2,i+1)=y0;
            result(3,i+1)=h*(i+1);
        end 

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
    f0=func(x0,y0);
        for i=1:n1  %快变量变化d次，慢变量再变化
            x1=x0+h1;
            for j=1:d
                y1=y0+h2*f0(2);
                f1=[x1,y1];
                f=IBMcalculate1(x0,y0, x1, f1,h2);
                f0=f;
                %x0=f(1);
                y0=f(2);            
                Y(end+1)=y0;
            end
            f=IBMcalculate2(x0,y0,x1,f,h1);
            x0=f(1);
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

function result = calculate(x0,y0, x1, f1, h)
    acc = 4;
    now = 0;
    f0=f1;
    iter=0;
    while now <= acc
        y1=f0(2);
        f1=func(x1,y1);
        f0=[x0+h*f1(1),y0+h*f1(2)];
        now =abs(log10(abs(f1(2)-f0(2))));
        iter=iter+1;
        if iter>5
            break;
        end
    end
    result = f0;
end

function result = IBMcalculate1(x0,y0, x1, f1,h2)
    acc = 4;
    now = 0;    
    f0=f1;   
    iter=0;
    while now <=acc
        y1=f0(2);
        f1=func(x1,y1);
        f0=[x0,y0+h2*f1(2);];     %更新y     
        now =abs(log10(abs(f1(2)-f0(2))));   
        iter=iter+1;
        if iter>5
            break;
        end
    end         
    result = f0;
end

function result = IBMcalculate2(x0,y0, x1, f1,h1)
    acc = 4;
    now = 0;    
    f0=f1;   
    iter=0;
    while now <=acc
        y1=f0(2);
        f1=func(x1,y1);
        f0=[x0+h1*f1(1),y0];   %更新x       
        now =abs(log10(abs(f1(1)-f0(1))));   
        iter=iter+1;
        if iter>5
            break;
        end
    end         
    result = f0;
end
