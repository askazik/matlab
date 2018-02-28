function test

clear
fileName = 'd:\!data\2018-02-13\20180216220227.frq';
D = classAmplitude(fileName);
disp(D)

h = figure('NumberTitle','off');
for j = 1:30000
    a = getUnit(D,j);
    frqs = D.frqs;
    heights = D.heights;
    set(h,'Name',int2str(a.gain));
    for i = 1:length(frqs)
        hh(i) = plot(abs(a.amplitudes(:,i)),heights);
        hold on
    end
    grid on
    pause(0.1)
    for i = 1:length(frqs)
        delete(hh(i))
    end
end

end