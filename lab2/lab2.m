disp("lab2");

function rez = var2(x, y)
  rez=x**2+1-3*y;
endfunction

function y = f(x)
  y = 1/27*(9.*x.**2-6*x+16*exp(-3.*x)+11);
endfunction

function [res_x,res_y] = RungeKut(func, step, n, x, y)
  res_x = [x];
  res_y = [y];
  for i =0:n
    k1=(feval(func,x,y));
    k2=feval(func,x+step/2,y+step/2*k1);
    k3=feval(func,x+step/2,y+step/2*k2);
    k4=feval(func,x+step,y+step*k3);
    y=y+step/6*(k1+2*k2+2*k3+k4);
    x=x+step;
    res_x = [res_x, x];
    res_y = [res_y, y];
  endfor
endfunction

function [res_x,res_y] = KutMerson(func, x0, X, y0, ep)
  res_x = [x0];
  res_y = [y0];
  x0=0;
  R=ep+1;
  x = x0;
  y = y0;
  h=0.1;
  do
    for i =0:h:X
      k1=h*feval(func,x,y);
      k2=h*feval(func,x+1/3*h,y+1/3*k1);
      k3=h*feval(func,x+1/3*h,y+1/6*k1+1/6*k2);
      k4=h*feval(func,x+1/2*h,y+1/8*k1+3/8*k3);
      k5=h*feval(func,x+h,y+1/2*k1-3/2*k3+2*k4);
      y1=y+1/2*k1-3/2*k3+2*k4;
      y=y+1/6*k1+2/3*k4+1/6*k5;
      R=0.2*abs(y1-y);
      if R>ep
        h=h/2;
        break
      endif
      x=x+h;
    res_x = [res_x, x];
    res_y = [res_y, y];
    endfor
  until (R>eps)
    
endfunction



x1 = linspace(0,1,101);
y1 = f(x1);



disp("RungeKut");
[x2,y2] = RungeKut('var2', 0.1, 9, 0, 1);
[x3,y3] = KutMerson('var2', 0,1,1,0.01);
subplot(1,2,1);
xlim([0,1]);
ylim([0.4,1]);
plot(x1,y1,'-b',x2,y2,'-g',x3,y3,'-m')

legend("Manual","RungeKut","KutMerson")

tspan = [0 1];
y0 = 1;
[x4,y4] = ode23(@(x,y) x**2+1-3*y, tspan, y0);
[x5,y5] = ode45(@(x,y) x**2+1-3*y, tspan, y0);
subplot(1,2,2)
plot(x1,y1,x4,y4,x5,y5)
legend("Manual","ode23","ode45")