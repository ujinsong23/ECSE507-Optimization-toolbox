function gradient = myGrad(f,x)
    gradient = zeros(size(x));
    h=1e-7;
    dimension = length(x);
    for d=1:dimension
        unit = zeros(size(x));
        unit(d)=1;
        gradient(d)=(f(x+h*unit)-f(x))/h;
    end
end
