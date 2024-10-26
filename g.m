function value=g(t) % 'val' means 'value'
global   x flag_BC flag_source
sp=200;
if flag_BC==1
    y=x(2:end-1);
else
    y=x;
end

val=exp(-sp*(y-0.5).^2).*10*(exp(-sp*(t-0.1).^2)+exp(-sp*(t-0.6).^2)+exp(-sp*(t-1.35).^2)+exp(-sp*(t-1.85).^2));
value=flag_source*val;
end
 