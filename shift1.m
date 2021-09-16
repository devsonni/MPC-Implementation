function [t0, x0, u0, xs] = shift1(T, t0, x0, u, f_u, xs)
st = x0;
con = u(1,:)';
f_value = f_u(st,con);
st = st+ (T*f_value);
x0 = full(st);

st_t = xs;
con_t = [15; 0.15];
f_t_value = [con_t(1)*cos(xs(3));con_t(1)*sin(xs(3));con_t(2)];
xs = xs + (T*f_t_value);
xs = full(xs);

t0 = t0 + T;
u0 = [u(2:size(u,1),:);u(size(u,1),:)];
end