function [x,xhat,fval] = penalty_barrier2d(type,f0,h1,h2,E,initial)
    stop=false;

    alpha_list = [10,20,50,100,200,500,1000,5000,10000];
    beta_list = [100,10,5,1,0.1,1e-2,1e-3,1e-4,1e-5];
    
    x(1,:)=initial;
    for j=1:length(alpha_list)
        alpha = alpha_list(j);
        beta = beta_list(j);
        if type=="penalty"
            if E
                F = @(x) f0(x(1),x(2))+alpha*(min(0,h1(x(1),x(2)))^2 + h2(x(1),x(2))^2);
            else
                F = @(x) f0(x(1),x(2))+alpha*(min(0,h1(x(1),x(2)))^2 + min(0,h2(x(1),x(2)))^2);
            end
        elseif type=="barrier"
            if E
                disp("barrier method doesn't allow equality constraint");
                stop=true;
                break;
            elseif h1(initial(1),initial(2))<0 | h2(initial(1),initial(2))<0
                disp("Initial point should be in feasible area to use barrier function");
                stop=true;
                break;
            else
                F = @(x) f0(x(1),x(2))+beta*(1/h1(x(1),x(2)) + 1/h2(x(1),x(2)));
            end
            
        else
            disp("type should be either 'penalty' or 'barrier");
            stop=true;
            break;
        end
        %[xhatstar,fval]=secant2d(F,x(j,:),false);
        [xhatstar,fval]=fminsearch(F,x(j,:),optimset('TolX',1e-10,'MaxFunEvals',10000,'MaxIter',10000));
        x(j+1,:)=xhatstar;
    end
    
    if stop==false
        xhat = x(j+1,:);
        fval = f0(x(j+1,1),x(j+1,2));
    else
        xhat = initial;
        fval = f0(initial(1),initial(2));
    end
    
end