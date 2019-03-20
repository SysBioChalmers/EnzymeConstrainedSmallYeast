
SAData = importdata('data/RxnAndSA.txt');
rxns = SAData.textdata;
specificActivity = SAData.data;
close all
histogramPostions = 1:0.25:4.5;
[numbers positions] = hist(log10(specificActivity), histogramPostions);  
hold on
baseColor = [0.18 0.33 0.52];
endColor = [0.32 0.77 1];


for i = 1:length(positions)
    colorValue = (positions(i)-min(positions))/max(positions);
    if colorValue < 0.5
        colorValue=colorValue/0.5;
    else
        colorValue = 1;
    end
    colorResult = baseColor + (endColor - baseColor) * colorValue;
    bar(positions(i), numbers(i), 0.2, 'FaceColor', colorResult);
end
hold off
xlabel('log10(specific activity)', 'FontSize',15,'FontName', 'Arial')
ylabel('Count', 'FontSize',15,'FontName', 'Arial')

set(gca,'FontSize',15,'FontName', 'Arial')

median(specificActivity)
mean(specificActivity)