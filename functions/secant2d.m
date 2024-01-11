function [x,xhat,fval] = secant2d(f, initial,plot,max_iter)
    interval = 0.05;
    if ~exist('max_iter','var')
        max_iter = 50000;
    end
    threshold = 0.0001;
   
    f2 = @(x) f(x(:,1),x(:,2));
    j=1;
    x(j,:) = initial;
    H(:,:,j) = eye(2);
    V(j) = f2(x(j,:));
    
    for j=1:max_iter
        
        stop=false;
        G(j,:) = myGrad(f2,x(j,:));
        V(j) = f2(x(j,:));
        
        
        if norm(G(j,:)) < threshold 
            disp(norm(G(j,:)));
            V(j) = f2(x(j,:));
            break
        end
        
        if j>1
            g(j-1,:)=G(j,:)-G(j-1,:);
            delx(j-1,:)=x(j,:)-x(j-1,:);
            
            numer = delx(j-1,:)' * delx(j-1,:);
            denom = delx(j-1,:) * g(j-1,:)';
            temp = H(:,:,j-1) + numer/denom;
            
            numer = (H(:,:,j-1)*g(j-1,:)')*(H(:,:,j-1)*g(j-1,:)')';
            denom = g(j-1,:)*H(:,:,j-1)*g(j-1,:)';
            H(:,:,j) = temp - numer/denom;
        end
        
        
        s(j,:)=-H(:,:,j)*G(j,:)';

        min_val = V(j);
        w(j)=0;
        for ww=interval*0.1*(-500:1000)
            if f2(x(j,:)+ww*s(j,:))<min_val
                min_val=f2(x(j,:)+ww*s(j,:));
                w(j)=ww;
            end
        end
        
        if w(j)==0
            stop=true;
            break
        end
        
        if j<max_iter
            x(j+1,:)=x(j,:)+w(j)*s(j,:);
        end
        
    end
    
    if stop 
        disp("w cannot be found");
        fprintf('number of iteration so far: %d \n',j);
        fprintf('current value : %.2f \n',V(j));
        disp(['current location: ', num2str(x(j,:))]);
    elseif j==max_iter
        disp("minimum could not be found within maximum iteration of our algorithm");
        fprintf('number of iteration so far: %d \n',j);
        fprintf('value after last iteration : %.2f \n',V(j));
        disp(['last location: ', num2str(x(j,:))]);
    else
        fprintf('number of total iteration: %d \n',j);
        fprintf('local minimum : %.2f \n',V(j));
        disp(['location of local minimum : ', num2str(x(j,:))]);
        
    end
    
    if plot==true
        xx =min(x(:,1))-interval*10:interval:max(x(:,1))+interval*10;
        yy =min(x(:,2))-interval*10:interval:max(x(:,2))+interval*10;
        surf(xx,yy,f(xx,yy'),'FaceAlpha',.5,'EdgeColor','interp'); hold on;
        scatter3(x(:,1),x(:,2),V,"red",'filled');hold on;
        plot3(x(:,1),x(:,2),V,"red");
        hold off;
    end

    xhat = x(j,:);
    fval = V(j);

end
