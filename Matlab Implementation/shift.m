function [t0, x0, u0, xs] = shift1(T, t0, x0, u, f_u, xs, sc)
st = x0;
con = u(1,:)';
f_value = f_u(st,con);
st = st+ (T*f_value);
x0 = full(st);

%updating target
v = 15;
con_t = [v, 0];  % Linear and angular velocity of target

if (sc-1 >= 300)
    con_t = [v, -(pi / 2) / 24];  % right turn
end
if (sc-1 >= 360)
    con_t = [v, 0];  % 50 ahead
end
if (sc-1 >= 410)
    con_t = [v, (pi / 2) / 24];  % left ahead
end
if (sc-1 >= 470)
    con_t = [v, 0];  % 100 ahead
end
if (sc-1 >= 570)
    con_t = [v, ((11 * pi) / 18) / 12];  % left 110 degree
end
if (sc-1 >= 630)
    con_t = [v, 0];  % 150 ahead
end
if (sc-1 >= 780)
    con_t = [v, ((7 * pi) / 18) / 12];  % 100 ahead
end
if (sc-1 >= 840)
    con_t = [v, 0];  % 100 ahead
end
if (sc-1 >= 940)
    con_t = [v, -(3 * pi / 18) / 12];  % right ahead
end
if (sc-1 >= 1000)
    con_t = [v, 0];  % 100 aheade
end
if (sc-1 >= 1100)
    con_t = [v, (3 * pi / 18) / 12];  % left ahead
end
if (sc-1 >= 1160)
    con_t = [v, 0];  % 100+77 ahead
end
if (sc-1 >= 1335)
    con_t = [v, (pi / 2) / 12];  % left ahead
end
if (sc-1 >= 1395)
    con_t = [v, 0];  % 100+35 ahead
end
if (sc-1 >= 1535)
    con_t = [v, (pi / 2) / 12];  % 100 ahead
end

f_t_value = [con_t(1)*cos(xs(3));con_t(1)*sin(xs(3));con_t(2)];
xs = xs + (T*f_t_value);
xs = full(xs);

t0 = t0 + T;
u0 = [u(2:size(u,1),:);u(size(u,1),:)];
end
