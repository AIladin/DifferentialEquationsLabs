disp("Lab4");

function rez = fun(x,t)
rez = t+x-1;
endfunction

function rez = solved(t)
rez = exp(t)-t;
endfunction

function [x, t] = miln(func, a, b, n, x0)
    h = (b - a) / n;
    x(1) = x0; 
    xp(1) = x(1);
    t = 1:n+1;
    for i = 1:n+1 
        t(i) = a + (i-1)*h;
    end
    
    for i = 2:4
        K1 = feval(func,t(i-1), x(i - 1));
        K2 = feval(func,t(i-1) + h/2, x(i-1) + h/2 * K1);
        K3 = feval(func,t(i-1) + h/2, x(i-1) + h/2 * K2);
        K4 = feval(func,t(i-1) + h, x(i-1) + h * K3);
        delta = h/6 * (K1 + 2*K2 + 2*K3 + K4);
        x(i) = x(i-1) + delta;
        xp(i) = x(i);
    end
    for i = 4:n 
        xp(i+1) = x(i-3)+4*h/3*(2*feval(func,t(i-2),x(i-2))-feval(func,t(i-1),x(i-1))+feval(func,t(i),x(i)));
        m = xp(i+1) + 28/29 * (x(i) - xp(i));
        x(i+1) = x(i-1)+h/3*(feval(func,t(i-1),x(i-1))+4*feval(func,t(i),x(i))+feval(func,t(i+1),m));
    end
end


[x1, t1] = miln("fun",0, 1, 100, 1);
t2 = 0:0.01:1;   
x2 = solved(t2);
plot(t1, x1, 'r', t2, x2, 'ck.')
legend("Miln's", "manual")

