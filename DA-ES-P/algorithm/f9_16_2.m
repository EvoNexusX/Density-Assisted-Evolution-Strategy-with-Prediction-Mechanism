max_run = 15;
for func = 9 : 16
   delete(gcp('nocreate'));
   parpool('local',max_run);
   spmd(max_run)
       disp(func),disp(labindex+15);
       DA_ES_P(func, labindex+15);
   end
   delete(gcp('nocreate'));
end
 