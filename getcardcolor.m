function cardcolor = getcardcolor(img, mask)

    % Extract the individual red, green, and blue color channels.
    redChannel = img(:, :, 1);
    greenChannel = img(:, :, 2);
    blueChannel = img(:, :, 3);
    % Get means
    meanR = mean(redChannel(mask));
    meanG = mean(greenChannel(mask));
    meanB = mean(blueChannel(mask));
    cardcolor = [meanR,meanG,meanB];