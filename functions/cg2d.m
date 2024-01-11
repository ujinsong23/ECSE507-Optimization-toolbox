function [x,xhat,fval] = cg2d(f, initial, plot,max_iter)
    interval = 0.1;
    if ~exist('max_iter','var')
        max_iter = 50000;
    end
    threshold = 0.0001;
   
    f2 = @(x) f(x(:,1),x(:,2));
    j=1;
    x(j,:) = initial;
    beta(j) = 0;
    V(j) = f2(x(j,:));
    
    for j=1:max_iter
        G(j,:) = myGrad(f2,x(j,:));
        V(j) = f2(x(j,:));
        
        if norm(G(j,:)) < threshold 
            disp(norm(G(j,:)));
            V(j) = f2(x(j,:));
            break
        end
        
        if j==1
            s(j,:)=-G(j,:);
        else
            beta(j,:)=(G(j,:)-G(j-1,:))*G(j,:)' / norm(G(j-1,:))^2 ;
            s(j,:)= - G(j,:) + beta(j,:)*s(j-1,:);
        end
        
        min_val = V(j);
        w(j)=0;
        for ww=interval*(1:1000)
            if f2(x(j,:)+ww*s(j,:))<min_val
                min_val=f2(x(j,:)+ww*s(j,:));
                w(j)=ww;
            end
        end
        
        if j<max_iter
            x(j+1,:) = x(j)+w(j)*s(j,:);
        end
        
        
    end
    
    if j==max_iter
        disp("minimum could not be found within maximum iteration of our algorithm")
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

