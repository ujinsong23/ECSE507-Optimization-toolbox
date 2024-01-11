function [x,xhat,fval] = sgd2d_armijo(f, initial,plot,max_iter)
    interval = 0.05;
    if ~exist('max_iter','var')
        max_iter = 10000;
    end
    threshold = 0.0001;
    
    f2 = @(x) f(x(:,1),x(:,2));
    j=1;
    x(j,:) = initial;
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
       
        s(j,:)=-G(j,:);
        
        
        delta=1.3;
        for p=0:20
            ww=delta^p;
            V_bar=V(j)+0.5*ww*G(j,:)*s(j,:)';
            if f2(x(j,:)+ww*s(j,:)) >= V_bar
                break;
            end
        end
        
        mu=0.85;
        for q=0:20
            ww=delta^p * mu^q;
            V_bar=V(j)+0.5*ww*G(j,:)*s(j,:)';
            if f2(x(j,:)+ww*s(j,:)) <= V_bar
                break;
            end
        end

        if p==10 | q==10
            stop=true;
            break;
        end
        
        w(j)=ww;
        
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
