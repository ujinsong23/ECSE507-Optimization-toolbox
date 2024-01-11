function [x,xhat,fval] = sgd_stdquad(a,b,C,initial)
    f2=@(x) a+b'*x+0.5*x'*C*x;
    interval = 0.05;
    max_iter = 50000;
    threshold = 0.0001;
    
    x = zeros(6,max_iter);
    j = 1;
    x(:,j) = initial;
    
    for j=1:max_iter
        stop=false;
        V(j) = f2(x(:,j));
        G(:,j) = myGrad(f2,x(:,j));
        
        if norm(G(:,j)) < threshold 
            disp(norm(G(:,j)));
            break
        end
       
        s(:,j)=-G(:,j);
        
        w(j)=rand(1)*interval;
        k=0;
        while true
            if f2(x(:,j)+w(j)*s(:,j)) < V(j)
                break
            elseif k<20
                w(j)=rand(1)*interval;
            else
                stop = true;
                break;
            end
            k=k+1;
        end
        
        if stop
            break;
        end
        
        if j<max_iter
            x(:,j+1)=x(:,j)+w(j)*s(:,j);
        end
    end
    
    if stop 
        disp("w cannot be found");
    else
        fprintf('number of total iteration: %d \n',j);
        fprintf('local minimum : %.2f \n',V(j));
        disp(['location of local minimum : ', num2str(x(:,j)')]);
    end
    
    xhat = x(:,j);
    fval = V(j);

end