#lab1 Korzhuk Andrew
display("lab1")

function res = var2(x, y)
  res=(y**2-y)/x; 
endfunction
  
function [res_x, res_y] = Euler(func, step, n, x, y)
  res_x = [x];
  res_y = [y];
  for i =0:n
    y=y+step*feval(func,x,y);
    x=x+step;
    res_x = [res_x, x];
    res_y = [res_y, y];
  endfor
  res = y;
endfunction

function [res_x,res_y] = RungeKut(func, step, n, alf, x, y)
  res_x = [x];
  res_y = [y];
  for i =0:n
    y=y+step*((1-alf)*feval(func,x,y)+alf*feval(func,x+step/(2*alf)*feval(func,x,y),y+step/(2*alf)));
    x=x+step;
    res_x = [res_x, x];
    res_y = [res_y, y];
  endfor
endfunction

function [res_x,res_y] = ModifiedEuler(fun, step, n, x, y)
  [res_x,res_y] = RungeKut(fun, step, n, 1, x, y);
endfunction

function [res_x,res_y] = Hoin(fun, step, n, x, y)
  [res_x,res_y] = RungeKut( fun, step, n, 1/2, x, y);
endfunction

function y = f(x)
  y = 1./(x.+1);
endfunction

x1 = linspace(1,2,11)
y1 = f(x)

disp("Euler")
[x2,y2] = Euler('var2', 0.1, 9, 1, 1/2)

disp("ModifiedEuler")
[x3,y3] = ModifiedEuler("var2", 0.1, 9, 1, 1/2)

disp("Hoin")
[x4,y4] = Hoin("var2", 0.1, 9, 1, 1/2)
plot(x1,y1,'-b',x2,y2,'-g',x3,y3,'-y',x4,y4,'-r')
legend("Manual","Euler","ModifiedEuler","Hoin")

