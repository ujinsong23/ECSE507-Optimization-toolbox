function [x,xhat,fval] = cg_stdquad(a,b,C,initial)
    f2=@(x) a+b'*x+0.5*x'*C*x;
    max_iter = 10000;
    threshold = 0.0001;
    
    x = zeros(6,max_iter);
    j = 1;
    x(:,j) = initial;
    
    for j=1:max_iter
        V(j) = f2(x(:,j));
        G(:,j) = myGrad(f2,x(:,j));
        
        if norm(G(:,j)) < threshold 
            disp(norm(G(:,j)));
            break
        end
       
        s(:,j) = -G(:,j);
        
        temp1 = -G(:,j)'*s(:,j);
        temp2 = s(:,j)'*C*s(:,j);
        w(j) = temp1/temp2;
        
        if j<max_iter
            x(:,j+1)=x(:,j)+w(j)*s(:,j);
        end
    end
    
    if j==max_iter 
        disp("w cannot be found");
    else
        fprintf('number of total iteration: %d \n',j);
        fprintf('local minimum : %.2f \n',V(j));
        disp(['location of local minimum : ', num2str(x(:,j)')]);
    end
    
    xhat = x(:,j);
    fval = V(j);

end