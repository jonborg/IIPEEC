K=[525 0 319.5;
    0 525 239.5;
    0 0 1];

depth1=imread('cardepth2.png');
im1=imread('car2.jpg');

depth2=imread('cardepth393.png');
im2=imread('car393.jpg');

xyz1=zeros(3, length(depth1(:)));
xyz2=zeros(3, length(depth2(:)));

u = repmat(1:640,480,1);
u = u(:)-K(1,3);
v = repmat((1:480)',640,1);
v=v(:)-K(2,3);
xyz1=zeros(length(u),3);
xyz2=zeros(length(u),3);


xyz1(:,3) = double(depth1(:))*0.001; % Convert to meters
xyz1(:,1) = (xyz1(:,3)/K(1,1)) .* u ;
xyz1(:,2) = (xyz1(:,3)/K(2,2)) .* v;
xyz2(:,3) = double(depth2(:))*0.001; % Convert to meters
xyz2(:,1) = (xyz2(:,3)/K(1,1)) .* u ;
xyz2(:,2) = (xyz2(:,3)/K(2,2)) .* v;

u1 =[329.1022  216.7007  324.5143 372.6864  352.0412 ]';
v1 =[111.8629  53.7925 80.2550 78.7848  194.9257 ]';
u2=[257.9910  224.7294  135.2670   42.3638  192.6147]';
v2=[87.6057  108.1876  124.3591  168.4632  130.9747]';

indice1=zeros(1,5);
indice2=zeros(1,5);
pontos1=zeros(5,3);
pontos2=zeros(5,3);

iterations=0;
erro=Inf;

for i=1:5
        indice1(i)=480*floor(u1(i))+floor(v1(i));
        indice2(i)=480*floor(u2(i))+floor(v2(i));

        pontos1(i,:)=xyz1(indice1(i),:);
        pontos2(i,:)=xyz2(indice2(i),:)';
end
R=eye(3);
T=zeros(3,1);

while (erro>0.1 & iterations<10)

    
    distancias=Inf(1,5);
    registo=[1 2 3 4 5; 0 0 0 0 0];
    taken=zeros(1,5);
    aux=zeros(5,3);
    for i=1:5
        for j=1:5
            if find(taken==j)
                continue;
            else
                if norm(pontos1(i,:)-pontos2(j,:))<distancias(i);
                    distancias(i)=norm(pontos1(i,:)-pontos2(j,:));
                    registo(2,i)=j;
                end
            end
        taken(i)=registo(2,i);
        end
        aux(i,:)=pontos2(registo(2,i),:);
    end
    pontos2=aux;
    
    centroid1=mean(pontos1);
    cdistancia1=pontos1-repmat(centroid1,5,1);
    centroid2=mean(pontos2);
    cdistancia2=pontos2-repmat(centroid2,5,1);

    Covariance=cdistancia1'*cdistancia2;

    [U,S,V] = svd(Covariance);

    Raux=U*V';
    Taux=centroid1'-Raux*centroid2';
    

    xyz2=Raux*xyz2'+repmat(Taux,1,480*640);
    xyz2=xyz2';
    
    for i=1:5
        pontos2(i,:)=xyz2(indice2(i),:)';
    end
    
    R=Raux*R;
    T=Raux*T+Taux;
    iterations=iterations+1;
    erro=sum(norm(pontos1-pontos2));
end

cl1=reshape(im1,480*640,3);
cl2=reshape(im2,480*640,3);
 
p1=pointCloud(xyz1,'Color',cl1);
p2=pointCloud(xyz2,'Color',cl2);

[R T]

figure(1);
showPointCloud(p1);
hold on
showPointCloud(p2);