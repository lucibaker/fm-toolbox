function [cross_up_idx, cross_down_idx] = flux_track(smtracks, Zc)
% return indices of smtracks where particles cross up or down through the level Zc

cross_up_idx = [0; smtracks(1:end-1,2) < Zc & smtracks(2:end,2) > Zc & ...
    smtracks(2:end,5) == smtracks(1:end-1,5)];

cross_down_idx = [0; smtracks(1:end-1,2) > Zc & smtracks(2:end,2) < Zc & ...
    smtracks(2:end,5) == smtracks(1:end-1,5)];