format long;
a=0;
b=3000;
h=0.0001;
x0=2; %方程初值
y0=0;

for i=1:2
    h=h*i;%变步长
    result=Euler(a,b,h,x0,y0);
    x=result(1,:);
    y=result(2,:);
    t=result(3,:);
    plot(t,x);
    hold on;%在同一张图绘制曲线
end
h=get(gca,'children');
set(h(1),'color','red'); %使第一条曲线变为红色


%% 需求解常微分方程
function fx=func(x,y)  
    fx=[y(1);
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

function result = calculate(x0,y0, x1, f1, h)

    acc = 10;
    now = 10;
    f0=f1;
    
    while now >= acc
        y1=f0(2);
        f1=func(x1,y1);
        f0=[x0+h*f1(1),y0+h*f1(2)];
        now = log10(abs(f1(2)-f0(2)));
    end
    result = f0;
end
