function [data,fits,layers,predict_fits] = GetFitness(data,pro,pop_I,layers,predict_fits)

    %% Evaluate or predict the fitness of the populations
    if pro.env+1<50
        fits_I = pro.GetFits(pop_I);
        fits = fits_I;
        data = [data,fits];
    else
        if pro.env+1 == 50
            fits_I = pro.GetFits(pop_I);
            fits = fits_I;
            data = [data,fits];
            mu = mean(data,2);
            sigma = std(data,0,2);
            normal_data = (data-mu)./sigma;
            options = trainingOptions('adam', 'MaxEpochs', 150, 'MiniBatchSize', 50);
    
            net = trainNetwork(normal_data(:,1:40)', normal_data(:,41:50)', layers, options);
            layers = net.Layers;
            predict_fits = predict(net,normal_data(:,11:50)')';
            predict_fits = predict_fits.*sigma+mu;
    
        else
%                   fits = predict_fits(:,pro.env-49);
                  fits_I = pro.GetFits(pop_I(1:2000,:));
                  p_fits = predict_fits(:,pro.env+1-50);
                  [~,index1] =  sort(fits_I);
                  [~,rank1]  =  sort(index1);
                  [~,index2] =  sort(p_fits(1:2000,:));
                  [~,rank2]  =  sort(index2);
                  sums = length(index2)*length(index2)/2;

                  total = sum(abs(rank1-rank2));
                  ac = 1-total/sums;
                  fprintf("准确率:%.1f%%",ac*100);
                  if ac>=0.60
                      disp("预测成功！");
                      fits = [fits_I;p_fits(2001:end,:)];
                  else
                      fits = [fits_I;pro.GetFits(pop_I(2001:end,:))];
                  end
%                   fits = pro.GetFits(pop_I);
                  data = [data,fits];
%                 fits_I = pro.GetFits(pop_I(1:2000,:));
%                 %               fits = fits_I;
%                 [results,nets,layers,mu,sigma] = net_train(nets,data,fits_I,pop_I,pro,layers,mu,sigma);
%                 %             results = fits_I';
%                 data = [data,results];
%                 fits = results;
        end
    end
 
end