max_run = 15;
for func = 21 : 24
   delete(gcp('nocreate'));
   parpool('local',max_run);
   spmd(max_run)
       disp(func),disp(labindex);
       DA_ES_P(func, labindex);
   end
   delete(gcp('nocreate'));
end
 