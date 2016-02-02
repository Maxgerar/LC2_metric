#include "Lc2.h"


//constructeur
Lc2::Lc2()
{
    /*****************************
     * SETTING DES IMAGES
     ****************************/
    cout<<"creation objet lc2"<<endl;
	//Clone l'image avant de la convertit en itk


	//Convertit l'image IRM de mitk en image itk
//	mitk::ImageToItk<ImageType>::Pointer toItkFilter =mitk::ImageToItk<ImageType>::New();
//    toItkFilter->SetInput(imageclone);
//    toItkFilter->Update();
//	 ImageType::Pointer itkImage = toItkFilter->GetOutput();
	
	 //transfert des geometries pour construire les matrices
//	geo_irm_=geo_irm;
//	geo_us_=geo_us;
//	image_irm_=image_irm;
//	image_us_=image_us;
	 
	 //Normalisation de matrice en fonction de l'image
//	Normalisationimage maximum;
//	maximum.calculemax(image_us,0);
//	maximum.calculemax(image_irm,1);
//	maximum.calculefacteur();
//	float facteurconversion=maximum.getfacteur();
//	cout<<"facteur intensité"<<endl;
//	cout<<facteurconversion<<endl;

	//calcule le gradient l'image

	//Reconvertit en image mitk
//	 mitk::Image::Pointer mitkgradient = mitk::ImportItkImage(gradientimage);

//	 //variable pour le calcul final de lc2
//	 double variancesum=0;
//	 double variancelc2poids=0;

    //initialisation des doubles
    variancesum =0;
    lc2varsum = 0;
    lc2final = 0;
	
}

void Lc2::setImages(ImageType::Pointer image_IRM, ImageType::Pointer image_US)
{
    
    
    //images IRM et US
    im_IRM = image_IRM;
    im_US = image_US;
    
    Gradient gradient;
    gradient.compute_gradient(image_IRM);
    im_grad_IRM = gradient.getImageGradient();
}

void Lc2::compute_param()
{
    cout<<"calcul de la lc2"<<endl;
    //setting du neighbourhood iterator
    // s=3
    ImageType::SizeType radius;
    radius[0] =3;
    radius[1] =3;
    radius[2] =3;
    
    //on itere sur toute l'image
    ImageType::RegionType accessibleImagePart;
    ImageType::RegionType::IndexType startIndex;
    startIndex[0] =3;
    startIndex[1] =3;
    startIndex[2] =3;
    ImageType::RegionType::SizeType sizeAccessible;
    ImageType::SizeType sizeIm = im_US->GetLargestPossibleRegion().GetSize();
    sizeAccessible[0] = sizeIm[0]-6;
    sizeAccessible[1] = sizeIm[1]-6;
    sizeAccessible[2] = sizeIm[2]-6;
    
    accessibleImagePart.SetIndex(startIndex);
    accessibleImagePart.SetSize(sizeAccessible);
    
    cout<<"index init de la region : "<<accessibleImagePart.GetIndex()<<endl;
    
    
    NIteratorType iterator(radius,im_US, accessibleImagePart); //im_US->GetLargestPossibleRegion()
    
    //placement au debut
    iterator.GoToBegin();
    // taille du voisinnage
    double m = (2*radius[0] + 1)^3;
    
    cout<<"index init iterator : "<<iterator.GetIndex(0)<<endl;
    
    
    //iteration sur l'image
    while(!iterator.IsAtEnd())
    {
        cout<<"new neighbourhood"<<endl;
        // on est donc sur un voisinnage diffŽrent ˆ chaque tour de la boucle while;
        //indice pour remplir les matrices
        int indice = 0;
        //statistiques sur le patch
        double mean=0;
        double variance=0;
        
       //matrices pour ce voisinnage
        VType U;
        //on rempli la matrice de 1, tout cela sera changŽ apres sauf la derniere colonne qui doit tre des 1
        MType M;
        M.Fill(1);
        
        //on construit un iterateur pour le neighbourhood
        ImageType::RegionType neighbourhood;
        ImageType::RegionType::IndexType start;
        //le debut de la region = le premier indice du masque
        start = iterator.GetIndex(0);
        cout<<"start index for neighbourhood iteration : "<<start<<endl;
        ImageType::RegionType::SizeType sizeN;
        sizeN[0] = 7;
        sizeN[1] = 7;
        sizeN[2] = 7;
        
        neighbourhood.SetIndex(start);
        neighbourhood.SetSize(sizeN);
        
        ImageConstIteratorType it(im_US,neighbourhood);
        
        it.GoToBegin();
        
        //parcours du voisinnage;
        
        while (!it.IsAtEnd()) {
            
            cout<<"parcours du voisinnage pour remplir les matrices"<<endl;
            //intensite US
            U[indice] = it.Get();
            
            //calcul de la moyenne d'intensitŽ sur le voisinnage
            mean = mean + it.Get();
            
            //on recupere le point dans l'espace correspondant pour aller le chercher dans l'IRM
            ImageType::IndexType indexUS = it.GetIndex();
            ImageType::PointType pt;
            
            im_US->TransformIndexToPhysicalPoint(indexUS, pt);
            ImageType::IndexType indexIRM;
            
            //si le point est dans l'image IRM
            if(im_IRM->TransformPhysicalPointToIndex(pt, indexIRM))
            {
                M(indice,0) = im_IRM->GetPixel(indexIRM);
                M(indice,1) = im_grad_IRM->GetPixel(indexIRM);
                
            }
            
            else
                //si on essaye de comparer a un element hors de l'image IRM
            {
                cout<<"en dehors de l'image IRM !"<<endl;
                //on met des 0 dans M !
                M(indice,0) = 0;
                M(indice,1) = 0;
            }
            
            indice++;
            ++it;
        }
        
        mean = mean/m;
        cout<<"moyenne : "<<mean<<endl;
        
        //calcul de la variance sur ce patch
        it.GoToBegin();
        cout<<"calcul de la variance pour patch"<<endl;
        while (!it.IsAtEnd()) {
            variance = variance + ((it.Get()- mean)*(it.Get() - mean));
            ++it;
        }
        
        variance = variance/m;
        cout<<"Variance : "<<variance<<endl;
        //on garde la valeur de variance pour ce voisinnage
        variances.push_back(variance);
        variancesum+= variance;
        
        cout<<"calcul des parametres par resol matricielle"<<endl;
        //matrice pour recuperer les params pour ce patchs = les coefficients de le relation lineaire entre US et IRM pour ce patch
        PType param;
        //MATRICES
        
        //affichage de M et U
        cout<<"matrice M : "<<endl;
        cout<<M<<endl;
        
        cout<<"Vecteur U : "<<endl;
        cout<<U<<endl;
       
        M3Type MTM = M3Type::Matrix(M.GetTranspose()*M.GetVnlMatrix());
        
        //inversion
        //c'est ici qu'il faut faire la gestion de det =0 !
        
        try {
            MTM = MTM.GetInverse();
            // c = (M^T*M)^-1 * M^T*U
            param = MTM*M.GetTranspose()*U;
        } catch (itk::ExceptionObject &e) {
            cerr<<"Matrix det is null"<<endl;
            cerr<<e<<endl;
            //ici c'est pcq ce sont des patchs unifomrme qu'il y a erreur
            param[0] = 1 ;
            param[1] = 0;
            param[2] = 0;
            
            
        }
        
        
        cout<<"parametres : "<<endl;
        cout<<"alpha : "<<param[0]<<endl;
        cout<<"beta : "<<param[1]<<endl;
        cout<<"gamma : "<<param[2]<<endl;
       
        //on enregistre cette valeur de params
        vector<double> paramvalue;
        paramvalue.clear();
        paramvalue.push_back(param[0]);
        paramvalue.push_back(param[1]);
        paramvalue.push_back(param[2]);
        
        parameters.push_back(paramvalue);
        
        
        //calcul de la LCI sur ce patch
        double sum = 0;
        
        cout<<"calcul LCI locale"<<endl;
        
        for (int j = 0; j<m;j++)
        {
            sum = sum +  ((U[j] - (param[0]*M(j,0) + param[1]*M(j,1) + param[2]))*(U[j] - (param[0]*M(j,0) + param[1]*M(j,1) + param[2])));
        }
        
        double lc2;
        if(variance != 0)
            lc2 = 1 - (sum/(m*variance));
        else
            lc2 = 1;
        
        cout<<"lc2 locale : "<<lc2<<endl;
        
        //ajout pondere par variance pour calcul de la LC2 totale
        lc2varsum = lc2varsum + (lc2*variance);
        
        LC2_i.push_back(lc2);
        
        
        //on incremente l'iterateur -> changer le voisinnage
        ++iterator;
        
    }
    
    //lc2 finale = moyenne ponderee
    
    lc2final = lc2varsum/variancesum;
    
//    //inscrire l'index du pixel central de la patch de départ sur l'US
//    itk::Index<3> idx = {{ 291,127,104}};
//    
//    
//    //section de l'image qui sera analysé à partir de l'index de départ (volume=longueur^3)
//    int longueur_section=3;
//    
//    for (int x=0;x<longueur_section;x++){
//        for (int y=0;y<longueur_section;y++){
//            for (int z=0;z<longueur_section;z++){
//                
//                //defini les pixels centrals de chaquep patch
//                itk::Index<3> idx_central;
//                idx_central[2]=idx[2]+x;
//                idx_central[1]=idx[1]+y;
//                idx_central[0]=idx[0]+z;
//                
//                //Le nombre de voisins du pixel central détermine la grandeur de la patch dans Matrice
//                int nombrevoisins=1;
//                
//                //construire matrice
//                Matricelc2 matricedelc2(nombrevoisins,idx_central,geo_irm_,geo_us_,image_irm_,image_us_,facteurconversion,mitkgradient);
//                //matricedelc2.affichemat();
//                //matricedelc2.affichevect();
//                
//                //construire matrice transverse
//                matricedelc2.domatricetransverse();
//                //matricedelc2.affichemattrans();
//                
//                //multiplication des matrices
//                matricedelc2.matricemultipli();
//                //matricedelc2.affichematavantinv();
//                
//                //inverse du résultats
//                matricedelc2.domatriceinverse();
//                //matricedelc2.affichematinv();
//                
//                matricedelc2.matricemultipli2();
//                //matricedelc2.affichematfinal();
//                
//                matricedelc2.calculecoefficient();
    
                
//                //affiche alpha beta gamma
//                matricedelc2.affichecoef();
//                
//                //calculer transformée de l'IRM
//                matricedelc2.calculevariance();
//                double variance=matricedelc2.getvariance();
//                matricedelc2.calculelc2local();
//                double lc2=matricedelc2.getlc2local();
//                //cout<<lc2<<endl;
//                
//                variancesum+=variance;
//                variancelc2poids+=variance*lc2;
//                
//            }
//        }
//    }
//    
//    lc2final=variancelc2poids/variancesum;
//    cout<<"lc2final"<<endl;
//    //cout<<lc2final<<endl;
    
}

double Lc2::getlc2(){
	return lc2final;
    
}

Lc2::MeasureType Lc2::GetValue(const TransformParametersType)
{
    compute_param();
    
    return lc2final;
    
}