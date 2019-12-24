th_1 = 0:1:50;
% x = 1:length(th_1);
% th_2 = 0:1:50;
% y = 1:length(th_1);
% th_3 = 0:1:50;
% z = 1:length(th_1);
lambdaVec = [1 3 5 7];
Q_gen=@(x,lambda) ((lambda.^x).*exp(-lambda)./factorial(x));
min_location = [0,0,0];
min=1;
for x=1:length(th_1)
    for y=x:length(th_1)
        for z=y:length(th_1)
            seg0 = 0:x; 
            seg1 = x:y;
            seg2 = y:z;
            seg3 = z:length(th_1);
            
            correct1 = Q_gen(seg0,lambdaVec(1));
            val1 = sum(correct1);
            correct2 = Q_gen(seg1,lambdaVec(2));
            val2 = sum(correct2);
            correct3 = Q_gen(seg2,lambdaVec(3));
            val3 = sum(correct3);
            correct4 = Q_gen(seg3,lambdaVec(4));
            val4 = sum(correct4);
            val = 1-0.25*(val1+val2+val3+val4);
            if (val<min)
                min = val;
                min_location=[x,y,z];
            end
        end
    end
end
