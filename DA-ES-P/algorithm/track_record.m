function [track] = track_reccord(groups)
    track = [];
    h =2;
    for i= 1:length(groups)
        track = [track;groups(i).OPTS.pop];
    end

end