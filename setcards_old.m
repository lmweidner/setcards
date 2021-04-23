
cardStats = struct();
r = [207,88,83];
g = [130,177,164];
p = [104,54,135];
rlab = [46,26];
glab = [-40,14];
plab = [37,-37];
colors = {'red','green','purple'};
%load image
img = imread('20210422_153623.jpg');

imflat = rgb2lab(reshape(img,size(img,1)*size(img,2),size(img,3)));
if max(imflat(:,1))<80
    img=imlocalbrighten(img,0.5);
    img=localcontrast(img,0.3,0.5);
end

%convert to grayscale and apply threshold
imggray = rgb2gray(img);
imgthresh = imggray>200;

%fill holes and extract boundaries
imgthresh2 = imfill(imgthresh,'holes');
figure(1)
imshow(imgthresh2)
[B,cardlab] = bwboundaries(imgthresh2,'noholes');
figure(2)
imshow(label2rgb(cardlab, @jet, [.5 .5 .5]))

c=0;
cards = {};
cardlabel = [];
for k = 1:length(B)
    
    boundsizes(k) = size(B{k},1);
    
end
boundthresh = 0.1*max(boundsizes);
for k =1:length(B)
    if size(B{k},1)>boundthresh %get regions that are larger than 50% of the largest region
        
        c=c+1;
        cardlabel(c) = k;
        cards{c} = B{k};
    end
end

%%

for i = 1:length(cards)
    vec=[];
    cardnum = i; %current card to analyze
    figure(3)
    mask = cardlab==cardlabel(cardnum); %mask of card in main image
    maskedRgbImage = bsxfun(@times, img, cast(mask, 'like', img));
    subplot(3,2,1)
    imshow(maskedRgbImage)
    
    %crop image to single card
    okind=find(mask>0);
    [ii,jj]=ind2sub(size(mask),okind);
    ymin=min(ii);ymax=max(ii);xmin=min(jj);xmax=max(jj);
    imCropped=imcrop(maskedRgbImage,[xmin,ymin,xmax-xmin+1,ymax-ymin+1]);
    subplot(3,2,2)
    imshow(imCropped)
    
    %kmeans segmentation of image
    %     [L1,Centers] = imsegkmeans(imCropped,3);
    %     B = labeloverlay(imCropped,L1);
    %     bwcrop = imbinarize(rgb2gray(imCropped));
    %     imshow(bwcrop)
    %     cc = bwconncomp(rgb2gray(imCropped))
    shapemask = createMask(imCropped);
    subplot(3,2,3)
    imshow(shapemask)
    title('Labeled Image')
    %     L4=L1(:);
    %     shapes2 = tabulate(L4(L4>0));
    %     [~,ind]= min(histcounts(L1(:))); %get cluster with the fewer number of pixels
    cardcolor = getcardcolor(imCropped,shapemask);
    if any(isnan(cardcolor))
        continue
    end
    %     for j = 1:max(L4)
    %         mask = L1==j; %mask of card shapes
    %         ctemp = getcardcolor(imCropped,mask);
    %         if sum(ctemp) > 10 && sum(ctemp) < 600
    %             cardcolor = ctemp;
    %         end
    %     end
    
    %display color of card
    cardcolorlab = rgb2lab(cardcolor./255);
    
    [~,idx]=min([pdist2(rlab,cardcolorlab(2:3),'euclidean'),pdist2(glab,cardcolorlab(2:3),'euclidean'),pdist2(plab,cardcolorlab(2:3),'euclidean')]);
    
    
    rawstats = regionprops(shapemask,'all');
    
    
    %     if length(rawstats)>3
    for ix = 1:length(rawstats)
        vec(ix) = rawstats(ix).Area;
    end
    threshold = 0.3*max(vec);
    numshapes = sum(vec>threshold);
    [~,largestshape] = max(vec);
    %     else
    
    %     end
    subplot(3,2,4)
    x = [0 1 1 0] ; y = [0 0 1 1] ;
    fill(x,y,cardcolor./255)
    title('card color')
    text(0.25,0.5,[num2str(numshapes),' shapes'],'fontsize',30)
    text(0.25,0.75,colors(idx),'fontsize',30)
    
    cardStats(i).num = numshapes;
    cardStats(i).color = colors(idx);
    
    %find number of shapes
    %     [B,L2] = bwboundaries(shapemask,'noholes');
    %     subplot(3,2,5)
    %     imshow(label2rgb(L2, @jet, [.5 .5 .5]))
    %     L3 = L2(:);
    %     shapes = tabulate(L3(L3>0));
    %     threshold = 0.5*max(shapes(:,2));
    %     numshapes = sum(shapes(:,2)>threshold);
    %     subplot(3,2,4)
    %     text(0.25,0.5,[num2str(numshapes),' shapes'],'fontsize',30)
    %
    % find shape and pattern
    %     [~,ind]=max(shapes(:,2)); %number of largest shape
    %     shapebound = B{ind};
    %     mask = L2 == ind;
    %     rawstats = regionprops(mask,'all');
    if length(rawstats)>1
        shape = rawstats(largestshape);
    else
        shape = rawstats;
    end
    
%     stat = (sum(sum(shape.FilledImage)) - sum(sum(shape.Image)))
    
    subplot(3,2,6)
    imshow(shape.FilledImage)
    stat = shape.Area/shape.FilledArea;
    if stat > 0.8
        subplot(3,2,4)
        text(0.5,0.25,['solid'],'fontsize',30)
        cardStats(i).inside = 'solid';
    elseif stat > 0.28 && stat < 0.8
        subplot(3,2,4)
        text(0.5,0.25,['hatched'],'fontsize',30)
        cardStats(i).inside = 'hatched';
    else
        subplot(3,2,4)
        text(0.5,0.25,['open'],'fontsize',30)
        cardStats(i).inside = 'open';
    end
    %     pause(1)
    
    %     shape.Perimeter/shape.Area
    %     shape.MajorAxisLength/shape.MinorAxisLength
    %     filledshape = regionprops(shape.FilledImage,'all');
    %     filledshape.ConvexArea - filledshape.Area
    
    subplot(3,2,5)
    H=[];
    T=[];
    R=[];
    lines=[];
    BW=[];
    % figure(101)
    % clf
    BW = edge(shape.FilledImage);
    [H,T,R] = hough(BW);
    peaks = houghpeaks(H,4,'threshold',0.3*max(H(:)));
    lines = houghlines(BW,T,R,peaks);
    max_len = 0;
    for k = 1:length(lines)
        xy = [lines(k).point1; lines(k).point2];
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
        hold on
        % Plot beginnings and ends of lines
        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    end
    hold off
    
    
    
    if (shape.ConvexArea - shape.FilledArea)/shape.ConvexArea > 0.1
        cardshape = 'squiggle';
        
        
        
        
        
    elseif length(lines) == 4
        line = squeeze(struct2cell(lines))';
        line2 = line(:,1:2);
        line3 = cell2mat(reshape(line2,8,1));
        linenums = [1,2,3,4,1,2,3,4];
        %     for ii = 1:8
        %         [k,dist] = dsearchn(line3,line3(ii,:));
        %
        %         if dist>0
        %             pair1(ii) = k;
        %         end
        %
        %     end
        d = squareform(pdist(line3));
        d = d + eye(8)*1e7;
        for ii = 1:8
            [~,pair1(ii)] = min(d(:,ii));
            pair2(ii) = linenums(pair1(ii));
        end
        pair3 = sort([linenums;pair2])';
        
        C = unique(pair3,'rows');
        C2 = tabulate(C(:));
        if sum(C2(:,2)==2)==4
            cardshape = 'diamond';
        else
            cardshape ='oval';
        end
    else
        cardshape = 'oval';
    end
    cardStats(i).shape = cardshape;
    
    cardStats(i).Centroid = rawstats.Centroid;
    cardStats(i).BB = rawstats.BoundingBox;
    

    subplot(3,2,4)
    text(0.25,0.95,cardshape,'fontsize',30)
    
    
    
    
    
    %crop image to single shape
    
    % okind=find(mask>0);
    % [ii,jj]=ind2sub(size(mask),okind);
    % ymin=min(ii);ymax=max(ii);xmin=min(jj);xmax=max(jj);
    % imCropped2=imcrop(imCropped,[xmin,ymin,xmax-xmin+1,ymax-ymin+1]);
    % figure(7)
    % imshow(imCropped2)
    
end

%%
figure(4)
clf
    imshow(img)
for ix = 1:length(cardStats)
    
   
    hold on
    boundary = cards{ix};
   plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2)
   text(mean(boundary(:,2)),mean(boundary(:,1)),...
       strcat(cardStats(ix).color,', ',cardStats(ix).inside,', ',cardStats(ix).shape,', ',num2str(cardStats(ix).num)),...
       'fontsize',12,'horizontalalignment','center')
    
end
hold off
