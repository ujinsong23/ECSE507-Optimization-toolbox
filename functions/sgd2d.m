function [x,xhat,fval] = sgd2d(f, initial,plot,max_iter)
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
        
        w(j)=rand(1)*interval;
        k=0;
        while true
            if f2(x(j,:)+w(j)*s(j,:)) < V(j)
                break
            elseif k<20
                w(j)=rand(1)*interval;
            elseif k<40
                w(j)=rand(1)*interval*0.1;
            elseif k<100
                w(j)=rand(1)*interval*20;
            else
                stop=true;
                break;
                
            end
            k=k+1;
        end
        
        if stop
            break;
        end
        
        if j<max_iter
            x(j+1,:)=x(j,:)+w(j)*s(j,:);
        end
    end
    
    if stop 
        disp("w cannot be found");
    elseif j==max_iter
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

