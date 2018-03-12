function a=am(v)       % Alpha for activation gating variable m in Na current
     E0m = v + 41;
     deltam = 1e-5; %mV
%      if (abs(E0m) < deltam)
%         a = 2000;
%      else
%         a = (200*E0m)/(1-exp(-0.1*E0m));
%      end

    %vectorize the command so can use in network
    ind_vec=abs(E0m)<deltam; %logical
    z_cont=(200*E0m)./(1-exp(-0.1*E0m)); %cont part
    z_cont(ind_vec)=1; %set this to any value; mult by 0 anyway
    a=2000*ind_vec + (1-ind_vec).*z_cont;
end


