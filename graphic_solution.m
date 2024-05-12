M = 5;
theta = 20;
    beta = theta+5:90;
    f = zeros(length(beta),1);
    for i = 1:length(beta)
        f(i) = 2*cotd(beta(i))*((M^2)*sind(beta(i))*sind(beta(i))-1)/((M^2)*(gemma+cosd(2*beta(i)))+2) - tand(theta);
    end
    plot(beta,f,'-r');
    hold on
theta = 30;
    beta = theta+1:90;
    f = zeros(length(beta),1);
    for i = 1:length(beta)
        f(i) = 2*cotd(beta(i))*((M^2)*sind(beta(i))*sind(beta(i))-1)/((M^2)*(gemma+cosd(2*beta(i)))+2) - tand(theta);
    end
    plot(beta,f,'-b');
    hold on
theta = 45;
    beta = theta+1:90;
    f = zeros(length(beta),1);
    for i = 1:length(beta)
        f(i) = 2*cotd(beta(i))*((M^2)*sind(beta(i))*sind(beta(i))-1)/((M^2)*(gemma+cosd(2*beta(i)))+2) - tand(theta);
    end
    plot(beta,f,'-c');
    hold on
    grid on
