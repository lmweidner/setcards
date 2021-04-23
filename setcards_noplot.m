%% process image and extract cards
brightnessparam = 0.8;
r = [207,88,83];
g = [130,177,164];
p = [104,54,135];
rlab = [46,26];
glab = [-40,14];
plab = [37,-37];
colors = {'red','green','purple'};
%load image
img = imread('20210422_153420.jpg');

imflat = rgb2lab(reshape(img,size(img,1)*size(img,2),size(img,3)));
if max(imflat(:,1))<80
    img=imlocalbrighten(img,0.5);
    img=localcontrast(img,0.3,0.5);
end

%convert to grayscale and apply threshold
imggray = rgb2gray(img);
imgthresh = imggray>255*brightnessparam;

%fill holes and extract boundaries
imgthresh2 = imfill(imgthresh,'holes');

[B,cardlab] = bwboundaries(imgthresh2,'noholes');


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

%% extract card properties
cardStats = {};
for i = 1:length(cards)
    
    vec=[];
    cardnum = i; %current card to analyze
    
    mask = cardlab==cardlabel(cardnum); %mask of card in main image
    maskedRgbImage = bsxfun(@times, img, cast(mask, 'like', img));

    
    %crop image to single card
    okind=find(mask>0);
    [ii,jj]=ind2sub(size(mask),okind);
    ymin=min(ii);ymax=max(ii);xmin=min(jj);xmax=max(jj);
    imCropped=imcrop(maskedRgbImage,[xmin,ymin,xmax-xmin+1,ymax-ymin+1]);
  
    
    shapemask = createMask(imCropped,brightnessparam);

    cardcolor = getcardcolor(imCropped,shapemask);
    if any(isnan(cardcolor))
        
        continue
    end
    
    %display color of card
    cardcolorlab = rgb2lab(cardcolor./255);
    
    [~,idx]=min([pdist2(rlab,cardcolorlab(2:3),'euclidean'),pdist2(glab,cardcolorlab(2:3),'euclidean'),pdist2(plab,cardcolorlab(2:3),'euclidean')]);
     cardStats{i,2} = colors{idx};
    
    rawstats = regionprops(shapemask,'all');
    
    
    %get number of shapes
    for ix = 1:length(rawstats)
        vec(ix) = rawstats(ix).Area;
    end
    threshold = 0.3*max(vec);
    numshapes = sum(vec>threshold);
    [~,largestshape] = max(vec);


    
    cardStats{i,1} = num2str(numshapes);
   
    

    if length(rawstats)>1
        shape = rawstats(largestshape);
    else
        shape = rawstats;
    end
    if shape.Area<10
        
        continue
    end
    
    %get shape fill
    stat = shape.Area/shape.FilledArea;
    if stat > 0.8

        cardStats{i,3} = 'solid';
    elseif stat > 0.28 && stat < 0.8

        cardStats{i,3} = 'hatched';
    else

        cardStats{i,3} = 'open';
    end

    %get shape type
    H=[];
    T=[];
    R=[];
    lines=[];
    BW=[];

    BW = edge(shape.FilledImage);
    [H,T,R] = hough(BW);
    peaks = houghpeaks(H,4,'threshold',0.3*max(H(:)));
    lines = houghlines(BW,T,R,peaks);
    max_len = 0;

    
    
    
    if (shape.ConvexArea - shape.FilledArea)/shape.ConvexArea > 0.1
        cardshape = 'squiggle';
        
        
        
        
        
    elseif length(lines) == 4
        line = squeeze(struct2cell(lines))';
        line2 = line(:,1:2);
        line3 = cell2mat(reshape(line2,8,1));
        linenums = [1,2,3,4,1,2,3,4];
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
    cardStats{i,4} = cardshape;
    
    cardStats{i,5} = rawstats.Centroid;
    cardStats{i,6} = rawstats.BoundingBox;
    
    
end
empt = any(cellfun(@isempty,cardStats),2);
cardStats(empt,:)=[];
cards(empt)=[];
%% plot extracted cards
figure(4)
clf
    imshow(img)
for ix = 1:size(cardStats,1)
    
   
    hold on
    boundary = cards{ix};
   plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2)
   text(mean(boundary(:,2)),mean(boundary(:,1)),...
       strcat(cardStats{ix,2},', ',cardStats{ix,3},', ',cardStats{ix,4},', ',cardStats{ix,1}),...
       'fontsize',12,'horizontalalignment','center')
    
end
hold off

%% find sets
matches = [];
nmatch = 0;
pnums = {'1','2','3'};
pcolors = {'red','green','purple'};
pinside = {'open','solid','hatched'};
pshape = {'diamond','oval','squiggle'};
for ix = 1:length(cardStats)
    oc = 1:length(cardStats);
    oc1 = oc(oc~=ix);
    for jx = oc1
        oc2 = oc1(oc1~=jx);
        card1 = cardStats(ix,1:4);
        card2 = cardStats(jx,1:4);
        card3={};
        same = [strcmp(card1{1},card2{1}),strcmp(card1{2},card2{2}),strcmp(card1{3},card2{3}),strcmp(card1{4},card2{4})];
        different = ~same;
        
        if different(1)
            num3 = pnums(~strcmp(pnums,card1{1}) & ~strcmp(pnums,card2{1}));
        else
            num3 = card1{1};
        end
        
        if different(2)
            color3 = pcolors(~strcmp(pcolors,card1{2}) & ~strcmp(pcolors,card2{2}));
        else
            color3 = card1{2};
        end
        
        if different(3)
            inside3 = pinside(~strcmp(pinside,card1{3}) & ~strcmp(pinside,card2{3}));
        else 
            inside3 = card1{3};
        end
        
        if different(4)
            shape3 = pshape(~strcmp(pshape,card1{4}) & ~strcmp(pshape,card2{4}));
        else
            shape3 = card1{4};
        end
        
        card3 = [num3,color3,inside3,shape3];
        
        
        
        for kx = oc2 %find if the third card is on the table
            
            match = length(intersect(cardStats(kx,1:4),card3))==4;
            if match %set found!
               nmatch = nmatch+1;
                matches(nmatch,:) = [ix,jx,kx];
                
            end
                
        end
    end
end
matches = unique(sort(matches,2),'rows');
if isempty(matches)
    disp('no sets found!')
end

%% show sets
figure(5)

for i = 1:size(matches,1)
    imshow(img)
    for j = 1:3
        
        hold on
        boundary = cards{matches(i,j)};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
        
    end
    hold off
    input('press enter to view next set')
end

