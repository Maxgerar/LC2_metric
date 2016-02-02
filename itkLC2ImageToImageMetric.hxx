/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkLC2ImageToImageMetric_hxx
#define itkLC2ImageToImageMetric_hxx

#include "itkLC2ImageToImageMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"


namespace itk
{
/**
 * Constructor
 */
template <typename TFixedImage, typename TMovingImage>
LC2ImageToImageMetric<TFixedImage, TMovingImage>
::LC2ImageToImageMetric()
{
  //m_SubtractMean = false;
    m_lc2final =0;
    variancesum=0;
    lc2varsum=0;
}
    
/*********************
 *Setting of images
 *********************/

template <typename TFixedImage, typename TMovingImage>
void LC2ImageToImageMetric<TFixedImage, TMovingImage>
::SetImages(typename TFixedImage::Pointer US, typename TMovingImage::Pointer IRM)
{
    this->m_MovingImage = IRM;
    this->m_FixedImage = US;
    
}
    

template<typename TFixedImage, typename TMovingImage>
    void LC2ImageToImageMetric<TFixedImage,TMovingImage>::SetFixed(typename TFixedImage::Pointer US)
    {
        this->m_FixedImage = US;
       
    }
    
template<typename TFixedImage, typename TMovingImage>
    void LC2ImageToImageMetric<TFixedImage,TMovingImage>::SetMoving(typename TMovingImage::Pointer IRM)
    {
        this->m_MovingImage =IRM;
    }
    
    


/*****
 * MASK FOR CROPPING
 *******/
    template <typename TFixedImage, typename TMovingImage>
    void
    LC2ImageToImageMetric<TFixedImage, TMovingImage>::ComputeMask()
    {
        cout<<"computing cropping mask"<<endl;
        //binarisation
        
        typename BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();
        thresholder->SetInput(this->m_FixedImage);
        thresholder->SetOutsideValue(255);
        thresholder->SetInsideValue(0);
        thresholder->SetLowerThreshold(0);
        thresholder->SetUpperThreshold(1);
        
        try {
            thresholder->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"error while binarizing US image"<<endl;
            cerr<<e<<endl;
        }
       
        MaskType::Pointer mask1 =thresholder->GetOutput() ;
      
        
        cout<<"done writing initial mask image"<<endl;
        
        cout<<"test closing"<<endl;
         //operation morphologiques
        

        kernelType::RadiusType radius;
        radius.Fill(3);
       
        kernelType kernel = kernelType::Ball(radius);
        cout<<"radius kernel : "<<kernel.GetRadius()<<endl;
        cout<<"kernel size : "<<kernel.GetSize()<<endl;
        
        CloserType::Pointer closer = CloserType::New();
        closer->SetInput(mask1);
        closer->SetKernel(kernel);
        closer->Update();
        
        m_mask = closer->GetOutput();
        
        //writing mask images
        
//        BinaryWriterType::Pointer writer3 = BinaryWriterType::New();
//        string out3 = "/Users/maximegerard/Documents/testmask.nii.gz";
//        itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//        writer3->SetInput(m_mask);
//        writer3->SetImageIO(io);
//        writer3->SetFileName(out3);
//        try {
//            writer3->Update();
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"error while writing image file"<<endl;
//            cerr<<e<<endl;
//            EXIT_FAILURE;
//        }
//        
        cout<<"done writing final mask image"<<endl;

        
    }
    
    
/** test the metric function
 */
    
template <typename TFixedImage, typename TMovingImage>
double
    LC2ImageToImageMetric<TFixedImage, TMovingImage>::test()
    {
        //resetting the varsum and lc2sum for new computation
        variancesum =0;
        lc2varsum = 0;
        
        //defining images
        FixedImageConstPointer fixedImage = this->m_FixedImage;
        
        if( !fixedImage )
        {
            itkExceptionMacro(<< "Fixed image has not been assigned");
        }

        
        
        //image qui bouge = image US
        MovingImageConstPointer movingImage = this->m_MovingImage;
        
        
        if(!movingImage)
        {
            itkExceptionMacro(<<"Moving image has not been assigned");
        }
        
        //gradient IRM
        
        cout<<"compute gradient of MRI image"<<endl;
        
        typename GradientFilterType::Pointer filterG = GradientFilterType::New();
        filterG->SetInput(movingImage);
        try {
            filterG->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"error while computing gradient image"<<endl;
            cerr<<e<<endl;
            EXIT_FAILURE;
        }
        
        typename TMovingImage::Pointer image_grad = filterG->GetOutput();
        
        //write image for test
        
//        typename WriterType::Pointer writer1 = WriterType::New();
//        string out1 = "/Users/maximegerard/Documents/testgrad.nii.gz";
//        NiftiImageIO::Pointer io = NiftiImageIO::New();
//        writer1->SetInput(image_grad);
//        writer1->SetImageIO(io);
//        writer1->SetFileName(out1);
//        try {
//            writer1->Update();
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"error while writing image file"<<endl;
//            cerr<<e<<endl;
//            EXIT_FAILURE;
//        }
//
//        cout<<"done writing gradient image"<<endl;
        
        ///////////////////////////////////////
        //CROPPING CONSIDERING MRI POSITION //
        //////////////////////////////////////
        
        
        //on s'interesse au voxel au milieu de chacun des face du volume MRI
        
        //get the limits of images in terms of indexes
        typename TMovingImage::SizeType sizeIRM = movingImage->GetLargestPossibleRegion().GetSize();
        cout<<"size IRM : "<<sizeIRM<<endl;
        typename TFixedImage::SizeType sizeUS = fixedImage->GetLargestPossibleRegion().GetSize();
        cout<<"size US : "<<sizeUS<<endl;
        
        typename TMovingImage::IndexType ind1;
        ind1[0] = 0;
        ind1[1] = sizeIRM[1]/2;
        ind1[2] = sizeIRM[2]/2;
        //cout<<"index IRM p1 : "<<ind1<<endl;
        typename TMovingImage::PointType p1;
        movingImage->TransformIndexToPhysicalPoint(ind1,p1);
        //cout<<"position spatiale p1 : "<<p1<<endl;;
        //eq US
        typename TFixedImage::IndexType indUS1;
        int indUS_infX;
        if(fixedImage->TransformPhysicalPointToIndex(p1,indUS1))
        {
            //cout<<"le point 1 est dans le volume US !"<<endl;
            //cout<<"indice US : "<<indUS1<<endl;
            indUS_infX = indUS1[0];
            
        }
        else
            indUS_infX = 0;
        
        typename TMovingImage::IndexType ind2;
        ind2[0] = sizeIRM[0]-1;
        ind2[1] = sizeIRM[1]/2;
        ind2[2] = sizeIRM[2]/2;
        //cout<<"indice irm p2 : "<<ind2<<endl;
        typename TMovingImage::PointType p2;
        movingImage->TransformIndexToPhysicalPoint(ind2,p2);
        //cout<<"position spatiale p2 : "<<p2<<endl;
        //eq US
        typename TFixedImage::IndexType indUS2;
        int indUS_supX;
        if(fixedImage->TransformPhysicalPointToIndex(p2,indUS2))
        {
            //cout<<"le point 2 est dans le volume US !"<<endl;
            //cout<<"indice US : "<<indUS2<<endl;
            indUS_supX = indUS2[0];
            
        }
        
        else
            indUS_supX = sizeUS[0]-1;

        
        
        typename TMovingImage::IndexType ind3;
        ind3[0] = sizeIRM[0]/2;
        ind3[1] = 0;
        ind3[2] = sizeIRM[2]/2;
        //cout<<"indice irm p3 : "<<ind3<<endl;
        typename TMovingImage::PointType p3;
        movingImage->TransformIndexToPhysicalPoint(ind3,p3);
        //cout<<"position spatiale p3 : "<<p3<<endl;
        //eq US
        typename TFixedImage::IndexType indUS3;
        int indUS_infY;
        if(fixedImage->TransformPhysicalPointToIndex(p3,indUS3))
        {
            //cout<<"le point 3 est dans le volume US !"<<endl;
            //cout<<"indice us : "<<indUS3<<endl;
            indUS_infY = indUS3[2];
            
        }
        else
            indUS_infY = 0;
        
        typename TMovingImage::IndexType ind4;
        ind4[0] = sizeIRM[0]/2;
        ind4[1] = sizeIRM[1]-1;
        ind4[2] = sizeIRM[2]/2;
        //cout<<"indice irm p4 : "<<ind4<<endl;
        typename TMovingImage::PointType p4;
        movingImage->TransformIndexToPhysicalPoint(ind4,p4);
        //cout<<"position spatiale p4 : "<<p4<<endl;
        //eq US
        typename TFixedImage::IndexType indUS4;
        int indUS_supY;
        if(fixedImage->TransformPhysicalPointToIndex(p4,indUS4))
        {
            //cout<<"le point 4 est dans le volume US !"<<endl;
            //cout<<"indice us : "<<indUS4<<endl;
            indUS_supY = indUS4[2];
            
        }
        else
            indUS_supY = sizeUS[2]-1;
        
        typename TMovingImage::IndexType ind5;
        ind5[0] = sizeIRM[0]/2;
        ind5[1] = sizeIRM[1]/2;
        ind5[2] = 0;
        //cout<<"indice irm p5 : "<<ind5<<endl;
        typename TMovingImage::PointType p5;
        movingImage->TransformIndexToPhysicalPoint(ind5,p5);
        //cout<<"position spatiale p5 : "<<p5<<endl;
        //eq US
        typename TFixedImage::IndexType indUS5;
        int indUS_infZ;
        if(fixedImage->TransformPhysicalPointToIndex(p5,indUS5))
        {
            //cout<<"le point 5 est dans le volume US !"<<endl;
            //cout<<"indice us : "<<indUS5<<endl;
            indUS_infZ = indUS5[1];
            
        }
        else
            indUS_infZ = 0;
        
        typename TMovingImage::IndexType ind6;
        ind6[0] = sizeIRM[0]/2;
        ind6[1] = sizeIRM[1]/2;
        ind6[2] = sizeIRM[2]-1;
        //cout<<"indice irm p6 : "<<ind6<<endl;
        typename TMovingImage::PointType p6;
        movingImage->TransformIndexToPhysicalPoint(ind6,p6);
        //cout<<"position spatiale p6 : "<<p6<<endl;
        //eq US
        typename TFixedImage::IndexType indUS6;
        int indUS_supZ;
        if(fixedImage->TransformPhysicalPointToIndex(p6,indUS6))
        {
            //cout<<"le point 6 est dans le volume US !"<<endl;
            indUS_supZ = indUS6[1];
            //cout<<"indice us : "<<indUS6<<endl;
            
        }
        else
            indUS_supZ = sizeUS[1]-1;
        
//        cout<<"test verification"<<endl;
//        cout<<"indice inf X "<<indUS_infX<<endl;
//        cout<<"indice sup X "<<indUS_supX<<endl;
//        cout<<"indice inf Y "<<indUS_infY<<endl;
//        cout<<"indice sup Y "<<indUS_supY<<endl;
//        cout<<"indice inf Z "<<indUS_infZ<<endl;
//        cout<<"indice sup Z "<<indUS_supZ<<endl;

        
        typename TFixedImage::IndexType startCropped;
        startCropped[0] = indUS_infX ;
        startCropped[1] = indUS_infZ;
        startCropped[2] = indUS_infY ;
        

        
        typename TFixedImage::IndexType endCropped;
        endCropped[0] = indUS_supX;
        endCropped[1] = indUS_supZ;
        endCropped[2] = indUS_supY;
        
        
        
        typename TFixedImage::SizeType sizeCropped;
        sizeCropped[0] = endCropped[0]-startCropped[0]+1;
        sizeCropped[1] = endCropped[1]-startCropped[1]+1;
        sizeCropped[2] = endCropped[2]-startCropped[2]+1;
        
        typename TFixedImage::RegionType regionCropped;
        regionCropped.SetIndex(startCropped);
        regionCropped.SetSize(sizeCropped);
        
        //Cropping de l'image originale
        
        typename ExtractorType::Pointer CroppingFilter = ExtractorType::New();
        CroppingFilter->SetExtractionRegion(regionCropped);
        CroppingFilter->SetInput(fixedImage);
        CroppingFilter->SetDirectionCollapseToIdentity();
        CroppingFilter->Update();
        typename TFixedImage::Pointer Fixed_Cropped = CroppingFilter->GetOutput();
        cout<<"verification image size : "<<Fixed_Cropped->GetLargestPossibleRegion().GetSize()<<endl;
        
//        //writing to verify
//        typename WriterType::Pointer writer4 = WriterType::New();
//        //itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//        string out4 = "/Users/maximegerard/Documents/testCroppingUS.nii.gz";
//        writer4->SetImageIO(io);
//        writer4->SetInput(Fixed_Cropped);
//        writer4->SetFileName(out4);
//        try {
//            writer4->Update();
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"error while writing croppped image"<<endl;
//            cerr<<e<<endl;
//        }
        
        //Cropping de l'image mask
        //pas necessaire on doit juste savoir a quelle region se limiter = region de cropped
        
//        typename BinaryExtractorType::Pointer BinaryCroppingFilter = BinaryExtractorType::New();
//        BinaryCroppingFilter->SetExtractionRegion(regionCropped);
//        BinaryCroppingFilter->SetInput(m_mask);
//        BinaryCroppingFilter->SetDirectionCollapseToIdentity();
//        BinaryCroppingFilter->Update();
//        MaskType::Pointer mask_Cropped = BinaryCroppingFilter->GetOutput();
        
//        //writing to verify
//        typename BinaryWriterType::Pointer writer5 = BinaryWriterType::New();
//        string out5 = "/Users/maximegerard/Documents/testCroppingMask.nii.gz";
//        writer5->SetImageIO(io);
//        writer5->SetInput(mask_Cropped);
//        writer5->SetFileName(out5);
//        try {
//            writer5->Update();
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"error while writing croppped mask image"<<endl;
//            cerr<<e<<endl;
//        }
    
        //on ne considere que les pixels non nuls de l'image US
        
        // On peut parcourir l'image avec le voisinnage mais ne faire les actions qui si le pixel eq dans le mask est a 1
       
        //////////////////////////////////////////
        // REGION ACCESSIBLE POUR NEIGHBORHOOD //
        ////////////////////////////////////////
        cout<<"defining accessible region of image"<<endl;
        //on itere sur toute l'image accessible -> celle pour laquelle le neighboorhood it ne va pas sortir de l'image
        typename TFixedImage::RegionType accessibleImagePart;
        typename TFixedImage::RegionType::IndexType startIndex = Fixed_Cropped->GetLargestPossibleRegion().GetIndex();
        startIndex[0] =startIndex[0]+3;
        startIndex[1] =startIndex[1]+3;
        startIndex[2] =startIndex[2]+3;
        typename TFixedImage::RegionType::SizeType sizeAccessible;
        typename TFixedImage::SizeType sizeIm = Fixed_Cropped->GetLargestPossibleRegion().GetSize();
        sizeAccessible[0] = sizeIm[0]-6;
        sizeAccessible[1] = sizeIm[1]-6;
        sizeAccessible[2] = sizeIm[2]-6;
        
        accessibleImagePart.SetIndex(startIndex);
        accessibleImagePart.SetSize(sizeAccessible);
        cout<<" start index region accessible : "<<startIndex<<endl;
        cout<<"taille region accessible : "<<accessibleImagePart.GetSize()<<endl;
        
        
        //TEST IMAGE ITERATOR RATHER THAN NITERATOR
        //CONST IT CAUSE WE DONT MODIFY IMAGE INTENSITIES
        
        
        //////////////////////////
        // ITERATIONS OVER IMAGE /
        /////////////////////////
        //to iterate over the accessible part of fixedImage
        ImageConstIteratorType US_it(fixedImage,accessibleImagePart);
        US_it.GoToBegin();
        

        cout<<"calcul de la lc2"<<endl;
        
        //to determine time of computation for lc2
        std::srand(time(NULL));
        std::time_t tbegin,tend;
        double texec = 0;
        tbegin = std::time(NULL);
        //size of neighbourhood
        int m = 343;
        
        //interpolateur pour image IRM et gradient
//        typename LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
//        interpolator->SetInputImage(movingImage);
//        
//        typename LinearInterpolatorType::Pointer gradInterpolator = LinearInterpolatorType::New();
//        gradInterpolator->SetInputImage(image_grad);
        
        
        while(!US_it.IsAtEnd())
        {
            //cout<<"US index under consideration : "<<US_it.GetIndex()<<endl;
            //on ne considere le voisinage que si le centre appartient à la region blanche du mask et s'il est a l'int de l'im IRM
            typename TMovingImage::PointType p;
            //get the equivalent in physical point to evaluate whether it's within the MRI image
            fixedImage->TransformIndexToPhysicalPoint(US_it.GetIndex(),p);
            //cout<<"Physical space coordinates of center of neighbourhood : "<<p<<endl;
            typename TMovingImage::IndexType i;
            //cout<<"verification mask -US : "<<int(m_mask->GetPixel(US_it.GetIndex()))<<endl;
            
            //we consider the neighbourhood only if it's center is in the actual US data and not outside of the MRI volume
            if(int(m_mask->GetPixel(US_it.GetIndex()))>0 && movingImage->TransformPhysicalPointToIndex(p,i))
            {
                //cout<<"neighbourhood in real US data"<<endl;
                // new neighborhood at each loop iteration
                
                //indice pour remplir les matrices
                int indice = 0;
                
                //MATRICES
                //matrices pour ce voisinnage
                VType U;
                //on rempli la matrice de 1, tout cela sera changé apres sauf la derniere colonne qui doit être des 1
                MType M;
                M.Fill(1);
                
                //STATS
                //statistiques sur le patch
                double mean=0;
                double variance=0;
                
                ///////////////////////////////////////////////////
                // ITERATION ON NEIGHBORHOOD : PARAM COMPUTATION
                /////////////////////////////////////////////////
                
                //define neighbourhood
                //here we define it as a region around the current pixel
                typename TFixedImage::RegionType neighbourhood;
                typename TFixedImage::RegionType::IndexType start;
                //le debut de la region = le premier indice du masque
                start = US_it.GetIndex();
                //le vrai debut est 3 pixel plus haut, plus a gauche et plus en profondeur
                start[0]= start[0]-3;
                start[1]= start[1]-3;
                start[2]= start[2]-3;
                
                //7-by-7 cube
                typename TFixedImage::RegionType::SizeType sizeN;
                sizeN[0] = 7;
                sizeN[1] = 7;
                sizeN[2] = 7;
                                
                neighbourhood.SetIndex(start);
                neighbourhood.SetSize(sizeN);
                
                ImageConstIteratorType it(fixedImage,neighbourhood);
                
                it.GoToBegin();
                
                //parcours du voisinnage;
                
                //NEIGHBORHOOD ITERATION
                while (!it.IsAtEnd())
                {
                
                    //intensite US
                    U[indice] = it.Get();
            
                    //calcul de la moyenne d'intensité US sur le voisinnage, mean computation is necessary to evaluate std dev
                    mean = mean + it.Get();
                
                    //on recupere le point dans l'espace correspondant pour aller le chercher dans l'IRM
                    typename TFixedImage::IndexType indexUS = it.GetIndex();
                    typename TFixedImage::PointType pt;
                
                    fixedImage->TransformIndexToPhysicalPoint(indexUS, pt);
                    //pt now contains the position in physica space of the considered voxel
                    
                    typename TMovingImage::IndexType indexIRM;
                
                    //si le point est dans l'image IRM
                    if(movingImage->TransformPhysicalPointToIndex(pt, indexIRM))
                    {
                        M(indice,0) = movingImage->GetPixel(indexIRM);
                        M(indice,1) = image_grad->GetPixel(indexIRM);
                        //test with linear interpolator
//                        M(indice,0) = interpolator->Evaluate(pt);
//                        M(indice,1) = gradInterpolator->Evaluate(pt);
                
                    }
                   
                
                    else//si on essaye de comparer a un element hors de l'image IRM
                    {
                        //cout<<"en dehors de l'image IRM !"<<endl;
                        //on met des 0 dans M !
                        // Ou on discard le voisinnage ?
                        M(indice,0) = 0;
                        M(indice,1) = 0;
                    }
                
                    indice++;
                    ++it;
                }
                
                mean = mean/m;
                //cout<<"moyenne : "<<mean<<endl;
                
                //STD DEV COMPUTATION
                //calcul de la variance sur ce patch
                it.GoToBegin();
                //cout<<"calcul de la variance pour patch"<<endl;
                while (!it.IsAtEnd()) {
                    variance = variance + ((it.Get()- mean)*(it.Get() - mean));
                    ++it;
                }
                
                variance = variance/m;
                //cout<<"Variance : "<<variance<<endl;
                //add to sum of variance for weighted average at the end to get global lc2
                variancesum+= variance;
                
                //MATRIX OPERATIONS
                //cout<<"calcul des parametres par resol matricielle"<<endl;
                //matrice pour recuperer les params pour ce patchs = les coefficients de le relation lineaire entre US et IRM pour ce patch
                PType param;
                //MATRICES
                
                //affichage de M et U
//                cout<<"matrice M : "<<endl;
//                cout<<M<<endl;
//                
//                cout<<"Vecteur U : "<<endl;
//                cout<<U<<endl;
                
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
                        cerr<<"spatial position of center of neighbourhood"<<p<<endl;
                        cerr<<"matric M"<<M<<endl;
                        cerr<<"metrice MTM"<<MTM<<endl;
                        //ici c'est pcq ce sont des patchs unifomrme qu'il y a erreur ?
                        param[0] = 1 ;
                        param[1] = 0;
                        param[2] = 0;
                
                
                    }
                
                                
                                
                //calcul de la LCI sur ce patch
                double sum = 0;
                
                ////////////
                //LOCAL LC2//
                ////////////
                //cout<<"calcul LCI locale"<<endl;
                
                //each pixel in neighbourhood has to be considered
                for (int j = 0; j<m;j++)
                {
                    sum = sum +  ((U[j] - (param[0]*M(j,0) + param[1]*M(j,1) + param[2]))*(U[j] - (param[0]*M(j,0) + param[1]*M(j,1) + param[2])));
                }
                                
                double lc2;
                if(variance != 0)
                {
                    lc2 = 1 - (sum/(m*variance));
                    //lc2 = sum/(m*variance);
                    
                }
                else
                {
                    cout<<"variance on patch is null"<<endl;
                    cout<<"neighbourhood voxel center : "<<p<<endl;
                    cout<<U<<endl;
                    lc2 = 0;
                    
                }
                
                                
                //cout<<"lc2 locale : "<<lc2<<endl;
                                
                //ajout pondere par variance pour calcul de la LC2 totale
                lc2varsum = lc2varsum + (lc2*variance);
                
                
                }
            
            
            ++US_it;

        }
        
        cout<<"done parcours image US"<<endl;
        tend = std::time(NULL);
        texec = std::difftime(tend,tbegin);
        
        cout<<"temps de parcours en s : "<<texec<<endl;
        
    //lc2 finale = moyenne ponderee
        
   
    m_lc2final = lc2varsum/variancesum;
        cout<<"lc2 globale : "<<m_lc2final<<endl;

        return m_lc2final;

        
    }


    
/**
 * Get the match Measure
 */
template <typename TFixedImage, typename TMovingImage>
typename LC2ImageToImageMetric<TFixedImage, TMovingImage>::MeasureType
LC2ImageToImageMetric<TFixedImage, TMovingImage>
::GetValue(const TransformParametersType & parameters) const
    //const TransformParametersType & parameters
{
    //setting des images
    
    //resetting the varsum and lc2sum for new computation
    double variancesum2 =0;
    double lc2varsum2 = 0;
    
    //defining images
    FixedImageConstPointer fixedImage = this->m_FixedImage;
    
    if( !fixedImage )
    {
        itkExceptionMacro(<< "Fixed image has not been assigned");
    }
    
    
    
    //image qui bouge = image US
    MovingImageConstPointer movingImage = this->m_MovingImage;
    
    
    if(!movingImage)
    {
        itkExceptionMacro(<<"Moving image has not been assigned");
    }
    
    //transformation with regard to TransformParameters
    
    //cout<<"test affichage param tsf : "<<parameters<<endl;
    TranslationTransformType::Pointer transform = TranslationTransformType::New();
    transform->SetParameters(parameters);
    cout<<"tsf parameters : "<<parameters<<endl;
    
    typename ResamplerType::Pointer resamplefilter = ResamplerType::New();
    resamplefilter->SetInput(movingImage);
    resamplefilter->SetSize(movingImage->GetLargestPossibleRegion().GetSize());
    resamplefilter->SetOutputSpacing(movingImage->GetSpacing());
    resamplefilter->SetOutputDirection(movingImage->GetDirection());
    //resamplefilter->SetTransform(transform->GetInverseTransform());
    resamplefilter->SetTransform(transform);

    cout<<"origin de l'image mobile before tsf : "<<movingImage->GetOrigin()<<endl;
    //typename TMovingImage::PointType origine = transform->TransformPoint(movingImage->GetOrigin());
    typename TMovingImage::PointType origine = (transform->GetInverseTransform())->TransformPoint(movingImage->GetOrigin());
    
    
//    origine[0] = origine[0]-parameters[0];
//    origine[1] = origine[1]-parameters[1];
//    origine[2] = origine[2]-parameters[2];
    
    cout<<"origin after tsf : "<<origine<<endl;
    
    resamplefilter->SetOutputOrigin(origine);
    try {
        resamplefilter->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while transforming moving image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    typename TMovingImage::Pointer movedImage = resamplefilter->GetOutput();

    //gradient IRM
    
    cout<<"compute gradient of MRI image"<<endl;
    
    typename GradientFilterType::Pointer filterG = GradientFilterType::New();
    filterG->SetInput(movedImage);
    try {
        filterG->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while computing gradient image"<<endl;
        cerr<<e<<endl;
        EXIT_FAILURE;
    }
    
    typename TMovingImage::Pointer image_grad = filterG->GetOutput();
    
    //write image for test
    
    //        typename WriterType::Pointer writer1 = WriterType::New();
    //        string out1 = "/Users/maximegerard/Documents/testgrad.nii.gz";
    //        NiftiImageIO::Pointer io = NiftiImageIO::New();
    //        writer1->SetInput(image_grad);
    //        writer1->SetImageIO(io);
    //        writer1->SetFileName(out1);
    //        try {
    //            writer1->Update();
    //        } catch (itk::ExceptionObject &e) {
    //            cerr<<"error while writing image file"<<endl;
    //            cerr<<e<<endl;
    //            EXIT_FAILURE;
    //        }
    //
    //        cout<<"done writing gradient image"<<endl;
    
    ///////////////////////////////////////
    //CROPPING CONSIDERING MRI POSITION //
    //////////////////////////////////////
    
    
    //on s'interesse au voxel au milieu de chacun des face du volume MRI
    
    //get the limits of images in terms of indexes
    typename TMovingImage::SizeType sizeIRM = movedImage->GetLargestPossibleRegion().GetSize();
    cout<<"size IRM : "<<sizeIRM<<endl;
    typename TFixedImage::SizeType sizeUS = fixedImage->GetLargestPossibleRegion().GetSize();
    cout<<"size US : "<<sizeUS<<endl;
    
    typename TMovingImage::IndexType ind1;
    ind1[0] = 0;
    ind1[1] = sizeIRM[1]/2;
    ind1[2] = sizeIRM[2]/2;
    //cout<<"index IRM p1 : "<<ind1<<endl;
    typename TMovingImage::PointType p1;
    movedImage->TransformIndexToPhysicalPoint(ind1,p1);
    //cout<<"position spatiale p1 : "<<p1<<endl;;
    //eq US
    typename TFixedImage::IndexType indUS1;
    int indUS_infX;
    if(fixedImage->TransformPhysicalPointToIndex(p1,indUS1))
    {
        //cout<<"le point 1 est dans le volume US !"<<endl;
        //cout<<"indice US : "<<indUS1<<endl;
        indUS_infX = indUS1[0];
        
    }
    else
        indUS_infX = 0;
    
    typename TMovingImage::IndexType ind2;
    ind2[0] = sizeIRM[0]-1;
    ind2[1] = sizeIRM[1]/2;
    ind2[2] = sizeIRM[2]/2;
    //cout<<"indice irm p2 : "<<ind2<<endl;
    typename TMovingImage::PointType p2;
    movedImage->TransformIndexToPhysicalPoint(ind2,p2);
    //cout<<"position spatiale p2 : "<<p2<<endl;
    //eq US
    typename TFixedImage::IndexType indUS2;
    int indUS_supX;
    if(fixedImage->TransformPhysicalPointToIndex(p2,indUS2))
    {
        //cout<<"le point 2 est dans le volume US !"<<endl;
        //cout<<"indice US : "<<indUS2<<endl;
        indUS_supX = indUS2[0];
        
    }
    
    else
        indUS_supX = sizeUS[0]-1;
    
    
    
    typename TMovingImage::IndexType ind3;
    ind3[0] = sizeIRM[0]/2;
    ind3[1] = 0;
    ind3[2] = sizeIRM[2]/2;
    //cout<<"indice irm p3 : "<<ind3<<endl;
    typename TMovingImage::PointType p3;
    movedImage->TransformIndexToPhysicalPoint(ind3,p3);
    //cout<<"position spatiale p3 : "<<p3<<endl;
    //eq US
    typename TFixedImage::IndexType indUS3;
    int indUS_infY;
    if(fixedImage->TransformPhysicalPointToIndex(p3,indUS3))
    {
        //cout<<"le point 3 est dans le volume US !"<<endl;
        //cout<<"indice us : "<<indUS3<<endl;
        indUS_infY = indUS3[2];
        
    }
    else
        indUS_infY = 0;
    
    typename TMovingImage::IndexType ind4;
    ind4[0] = sizeIRM[0]/2;
    ind4[1] = sizeIRM[1]-1;
    ind4[2] = sizeIRM[2]/2;
    //cout<<"indice irm p4 : "<<ind4<<endl;
    typename TMovingImage::PointType p4;
    movedImage->TransformIndexToPhysicalPoint(ind4,p4);
    //cout<<"position spatiale p4 : "<<p4<<endl;
    //eq US
    typename TFixedImage::IndexType indUS4;
    int indUS_supY;
    if(fixedImage->TransformPhysicalPointToIndex(p4,indUS4))
    {
        //cout<<"le point 4 est dans le volume US !"<<endl;
        //cout<<"indice us : "<<indUS4<<endl;
        indUS_supY = indUS4[2];
        
    }
    else
        indUS_supY = sizeUS[2]-1;
    
    typename TMovingImage::IndexType ind5;
    ind5[0] = sizeIRM[0]/2;
    ind5[1] = sizeIRM[1]/2;
    ind5[2] = 0;
    //cout<<"indice irm p5 : "<<ind5<<endl;
    typename TMovingImage::PointType p5;
    movedImage->TransformIndexToPhysicalPoint(ind5,p5);
    //cout<<"position spatiale p5 : "<<p5<<endl;
    //eq US
    typename TFixedImage::IndexType indUS5;
    int indUS_infZ;
    if(fixedImage->TransformPhysicalPointToIndex(p5,indUS5))
    {
        //cout<<"le point 5 est dans le volume US !"<<endl;
        //cout<<"indice us : "<<indUS5<<endl;
        indUS_infZ = indUS5[1];
        
    }
    else
        indUS_infZ = 0;
    
    typename TMovingImage::IndexType ind6;
    ind6[0] = sizeIRM[0]/2;
    ind6[1] = sizeIRM[1]/2;
    ind6[2] = sizeIRM[2]-1;
    //cout<<"indice irm p6 : "<<ind6<<endl;
    typename TMovingImage::PointType p6;
    movedImage->TransformIndexToPhysicalPoint(ind6,p6);
    //cout<<"position spatiale p6 : "<<p6<<endl;
    //eq US
    typename TFixedImage::IndexType indUS6;
    int indUS_supZ;
    if(fixedImage->TransformPhysicalPointToIndex(p6,indUS6))
    {
        //cout<<"le point 6 est dans le volume US !"<<endl;
        indUS_supZ = indUS6[1];
        //cout<<"indice us : "<<indUS6<<endl;
        
    }
    else
        indUS_supZ = sizeUS[1]-1;
    
    //        cout<<"test verification"<<endl;
    //        cout<<"indice inf X "<<indUS_infX<<endl;
    //        cout<<"indice sup X "<<indUS_supX<<endl;
    //        cout<<"indice inf Y "<<indUS_infY<<endl;
    //        cout<<"indice sup Y "<<indUS_supY<<endl;
    //        cout<<"indice inf Z "<<indUS_infZ<<endl;
    //        cout<<"indice sup Z "<<indUS_supZ<<endl;
    
    
    typename TFixedImage::IndexType startCropped;
    startCropped[0] = indUS_infX ;
    startCropped[1] = indUS_infZ;
    startCropped[2] = indUS_infY ;
    
    
    
    typename TFixedImage::IndexType endCropped;
    endCropped[0] = indUS_supX;
    endCropped[1] = indUS_supZ;
    endCropped[2] = indUS_supY;
    
    
    
    typename TFixedImage::SizeType sizeCropped;
    sizeCropped[0] = endCropped[0]-startCropped[0]+1;
    sizeCropped[1] = endCropped[1]-startCropped[1]+1;
    sizeCropped[2] = endCropped[2]-startCropped[2]+1;
    
    typename TFixedImage::RegionType regionCropped;
    regionCropped.SetIndex(startCropped);
    regionCropped.SetSize(sizeCropped);
    
    //Cropping de l'image originale
    
    typename ExtractorType::Pointer CroppingFilter = ExtractorType::New();
    CroppingFilter->SetExtractionRegion(regionCropped);
    CroppingFilter->SetInput(fixedImage);
    CroppingFilter->SetDirectionCollapseToIdentity();
    CroppingFilter->Update();
    typename TFixedImage::Pointer Fixed_Cropped = CroppingFilter->GetOutput();
    cout<<"verification image size : "<<Fixed_Cropped->GetLargestPossibleRegion().GetSize()<<endl;
    
    //        //writing to verify
    //        typename WriterType::Pointer writer4 = WriterType::New();
    //        //itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
    //        string out4 = "/Users/maximegerard/Documents/testCroppingUS.nii.gz";
    //        writer4->SetImageIO(io);
    //        writer4->SetInput(Fixed_Cropped);
    //        writer4->SetFileName(out4);
    //        try {
    //            writer4->Update();
    //        } catch (itk::ExceptionObject &e) {
    //            cerr<<"error while writing croppped image"<<endl;
    //            cerr<<e<<endl;
    //        }
    
    //Cropping de l'image mask
    //pas necessaire on doit juste savoir a quelle region se limiter = region de cropped
    
    //        typename BinaryExtractorType::Pointer BinaryCroppingFilter = BinaryExtractorType::New();
    //        BinaryCroppingFilter->SetExtractionRegion(regionCropped);
    //        BinaryCroppingFilter->SetInput(m_mask);
    //        BinaryCroppingFilter->SetDirectionCollapseToIdentity();
    //        BinaryCroppingFilter->Update();
    //        MaskType::Pointer mask_Cropped = BinaryCroppingFilter->GetOutput();
    
    //        //writing to verify
    //        typename BinaryWriterType::Pointer writer5 = BinaryWriterType::New();
    //        string out5 = "/Users/maximegerard/Documents/testCroppingMask.nii.gz";
    //        writer5->SetImageIO(io);
    //        writer5->SetInput(mask_Cropped);
    //        writer5->SetFileName(out5);
    //        try {
    //            writer5->Update();
    //        } catch (itk::ExceptionObject &e) {
    //            cerr<<"error while writing croppped mask image"<<endl;
    //            cerr<<e<<endl;
    //        }
    
    //on ne considere que les pixels non nuls de l'image US
    
    // On peut parcourir l'image avec le voisinnage mais ne faire les actions qui si le pixel eq dans le mask est a 1
    
    //////////////////////////////////////////
    // REGION ACCESSIBLE POUR NEIGHBORHOOD //
    ////////////////////////////////////////
    cout<<"defining accessible region of image"<<endl;
    //on itere sur toute l'image accessible -> celle pour laquelle le neighboorhood it ne va pas sortir de l'image
    typename TFixedImage::RegionType accessibleImagePart;
    typename TFixedImage::RegionType::IndexType startIndex = Fixed_Cropped->GetLargestPossibleRegion().GetIndex();
    startIndex[0] =startIndex[0]+3;
    startIndex[1] =startIndex[1]+3;
    startIndex[2] =startIndex[2]+3;
    typename TFixedImage::RegionType::SizeType sizeAccessible;
    typename TFixedImage::SizeType sizeIm = Fixed_Cropped->GetLargestPossibleRegion().GetSize();
    sizeAccessible[0] = sizeIm[0]-6;
    sizeAccessible[1] = sizeIm[1]-6;
    sizeAccessible[2] = sizeIm[2]-6;
    
    accessibleImagePart.SetIndex(startIndex);
    accessibleImagePart.SetSize(sizeAccessible);
    cout<<" start index region accessible : "<<startIndex<<endl;
    cout<<"taille region accessible : "<<accessibleImagePart.GetSize()<<endl;
    
    
    //TEST IMAGE ITERATOR RATHER THAN NITERATOR
    //CONST IT CAUSE WE DONT MODIFY IMAGE INTENSITIES
    
    
    //////////////////////////
    // ITERATIONS OVER IMAGE /
    /////////////////////////
    //to iterate over the accessible part of fixedImage
    ImageConstIteratorType US_it(fixedImage,accessibleImagePart);
    US_it.GoToBegin();
    
    
    cout<<"calcul de la lc2"<<endl;
    
    //to determine time of computation for lc2
    std::srand(time(NULL));
    std::time_t tbegin,tend;
    double texec = 0;
    tbegin = std::time(NULL);
    //size of neighbourhood
    int m = 343;
    
    //interpolateur pour image IRM et gradient
    //        typename LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
    //        interpolator->SetInputImage(movedImage);
    //
    //        typename LinearInterpolatorType::Pointer gradInterpolator = LinearInterpolatorType::New();
    //        gradInterpolator->SetInputImage(image_grad);
    
    
    while(!US_it.IsAtEnd())
    {
        //cout<<"US index under consideration : "<<US_it.GetIndex()<<endl;
        //on ne considere le voisinage que si le centre appartient à la region blanche du mask et s'il est a l'int de l'im IRM
        typename TMovingImage::PointType p;
        //get the equivalent in physical point to evaluate whether it's within the MRI image
        fixedImage->TransformIndexToPhysicalPoint(US_it.GetIndex(),p);
        //cout<<"Physical space coordinates of center of neighbourhood : "<<p<<endl;
        typename TMovingImage::IndexType i;
        //cout<<"verification mask -US : "<<int(m_mask->GetPixel(US_it.GetIndex()))<<endl;
        
        //we consider the neighbourhood only if it's center is in the actual US data and not outside of the MRI volume
        if(int(m_mask->GetPixel(US_it.GetIndex()))>0 && movedImage->TransformPhysicalPointToIndex(p,i))
        {
            //cout<<"neighbourhood in real US data"<<endl;
            // new neighborhood at each loop iteration
            
            //indice pour remplir les matrices
            int indice = 0;
            
            //MATRICES
            //matrices pour ce voisinnage
            VType U;
            //on rempli la matrice de 1, tout cela sera changé apres sauf la derniere colonne qui doit être des 1
            MType M;
            M.Fill(1);
            
            //STATS
            //statistiques sur le patch
            double mean=0;
            double variance=0;
            
            ///////////////////////////////////////////////////
            // ITERATION ON NEIGHBORHOOD : PARAM COMPUTATION
            /////////////////////////////////////////////////
            
            //define neighbourhood
            //here we define it as a region around the current pixel
            typename TFixedImage::RegionType neighbourhood;
            typename TFixedImage::RegionType::IndexType start;
            //le debut de la region = le premier indice du masque
            start = US_it.GetIndex();
            //le vrai debut est 3 pixel plus haut, plus a gauche et plus en profondeur
            start[0]= start[0]-3;
            start[1]= start[1]-3;
            start[2]= start[2]-3;
            
            //7-by-7 cube
            typename TFixedImage::RegionType::SizeType sizeN;
            sizeN[0] = 7;
            sizeN[1] = 7;
            sizeN[2] = 7;
            
            neighbourhood.SetIndex(start);
            neighbourhood.SetSize(sizeN);
            
            ImageConstIteratorType it(fixedImage,neighbourhood);
            
            it.GoToBegin();
            
            //parcours du voisinnage;
            
            //NEIGHBORHOOD ITERATION
            while (!it.IsAtEnd())
            {
                
                //intensite US
                U[indice] = it.Get();
                
                //calcul de la moyenne d'intensité US sur le voisinnage, mean computation is necessary to evaluate std dev
                mean = mean + it.Get();
                
                //on recupere le point dans l'espace correspondant pour aller le chercher dans l'IRM
                typename TFixedImage::IndexType indexUS = it.GetIndex();
                typename TFixedImage::PointType pt;
                
                fixedImage->TransformIndexToPhysicalPoint(indexUS, pt);
                //pt now contains the position in physica space of the considered voxel
                
                typename TMovingImage::IndexType indexIRM;
                
                //si le point est dans l'image IRM
                if(movedImage->TransformPhysicalPointToIndex(pt, indexIRM))
                {
                    M(indice,0) = movedImage->GetPixel(indexIRM);
                    M(indice,1) = image_grad->GetPixel(indexIRM);
                    //test with linear interpolator
                    //                        M(indice,0) = interpolator->Evaluate(pt);
                    //                        M(indice,1) = gradInterpolator->Evaluate(pt);
                    
                }
                
                
                else//si on essaye de comparer a un element hors de l'image IRM
                {
                    //cout<<"en dehors de l'image IRM !"<<endl;
                    //on met des 0 dans M !
                    // Ou on discard le voisinnage ?
                    M(indice,0) = 0;
                    M(indice,1) = 0;
                }
                
                indice++;
                ++it;
            }
            
            mean = mean/m;
            //cout<<"moyenne : "<<mean<<endl;
            
            //STD DEV COMPUTATION
            //calcul de la variance sur ce patch
            it.GoToBegin();
            //cout<<"calcul de la variance pour patch"<<endl;
            while (!it.IsAtEnd()) {
                variance = variance + ((it.Get()- mean)*(it.Get() - mean));
                ++it;
            }
            
            variance = variance/m;
            //cout<<"Variance : "<<variance<<endl;
            //add to sum of variance for weighted average at the end to get global lc2
            variancesum2+= variance;
            
            //MATRIX OPERATIONS
            //cout<<"calcul des parametres par resol matricielle"<<endl;
            //matrice pour recuperer les params pour ce patchs = les coefficients de le relation lineaire entre US et IRM pour ce patch
            PType param;
            //MATRICES
            
            //affichage de M et U
            //                cout<<"matrice M : "<<endl;
            //                cout<<M<<endl;
            //
            //                cout<<"Vecteur U : "<<endl;
            //                cout<<U<<endl;
            
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
                cerr<<"spatial position of center of neighbourhood"<<p<<endl;
                cerr<<"matric M"<<M<<endl;
                cerr<<"metrice MTM"<<MTM<<endl;
                //ici c'est pcq ce sont des patchs unifomrme qu'il y a erreur ?
                param[0] = 1 ;
                param[1] = 0;
                param[2] = 0;
                
                
            }
            
            
            
            //calcul de la LCI sur ce patch
            double sum = 0;
            
            ////////////
            //LOCAL LC2//
            ////////////
            //cout<<"calcul LCI locale"<<endl;
            
            //each pixel in neighbourhood has to be considered
            for (int j = 0; j<m;j++)
            {
                sum = sum +  ((U[j] - (param[0]*M(j,0) + param[1]*M(j,1) + param[2]))*(U[j] - (param[0]*M(j,0) + param[1]*M(j,1) + param[2])));
            }
            
            double lc2;
            if(variance != 0)
            {
                lc2 = 1 - (sum/(m*variance));
                //lc2 = sum/(m*variance);
                
            }
            else
            {
                cout<<"variance on patch is null"<<endl;
                cout<<"neighbourhood voxel center : "<<p<<endl;
                cout<<U<<endl;
                lc2 = 0;
                
            }
            
            
            //cout<<"lc2 locale : "<<lc2<<endl;
            
            //ajout pondere par variance pour calcul de la LC2 totale
            lc2varsum2 = lc2varsum2 + (lc2*variance);
            
            
        }
        
        
        ++US_it;
        
    }
    
    cout<<"done parcours image US"<<endl;
    tend = std::time(NULL);
    texec = std::difftime(tend,tbegin);
    
    cout<<"temps de parcours en s : "<<texec<<endl;
    
    //lc2 finale = moyenne ponderee
    
    
    double lc2final = lc2varsum2/variancesum2;
    cout<<"lc2 globale : "<<lc2final<<endl;
    

 
    

    MeasureType measure = lc2final;

  return measure;
}

/**
 * Get the Derivative Measure
 */
template <typename TFixedImage, typename TMovingImage>
void
LC2ImageToImageMetric<TFixedImage, TMovingImage>
::GetDerivative(const TransformParametersType & parameters,
                DerivativeType & derivative) const
{
  
}

/*
 * Get both the match Measure and theDerivative Measure
 */
template <typename TFixedImage, typename TMovingImage>
void
LC2ImageToImageMetric<TFixedImage, TMovingImage>
::GetValueAndDerivative(const TransformParametersType & parameters,
                        MeasureType & value, DerivativeType  & derivative) const
{
    
}

template <typename TFixedImage, typename TMovingImage>
void
LC2ImageToImageMetric<TFixedImage, TMovingImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "LC2final: " << m_lc2final << std::endl;
}

} // end namespace itk

template <typename TFixedImage, typename TMovingImage>
void
itk::LC2ImageToImageMetric<TFixedImage, TMovingImage>
::operator=(const Self &)
{
    
}

template <typename TFixedImage, typename TMovingImage>
itk::LC2ImageToImageMetric<TFixedImage, TMovingImage>
::LC2ImageToImageMetric(const Self &)
{
    
}



#endif
