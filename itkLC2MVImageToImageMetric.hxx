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
#ifndef itkLC2MVImageToImageMetric_hxx
#define itkLC2MVImageToImageMetric_hxx

#include "itkLC2MVImageToImageMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkArray.h"
#include <vector>

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedImage, typename TMovingImage>
LC2MVImageToImageMetric<TFixedImage, TMovingImage>
::LC2MVImageToImageMetric()
{
    m_numberofValues = 10000;
    //m_Transform = EulerTransformType::New();
  
}
    
template <typename TFixedImage, typename TMovingImage>
void
    LC2MVImageToImageMetric<TFixedImage, TMovingImage>::SetImages(typename TFixedImage::Pointer US, typename TMovingImage::Pointer IRM)
{
    m_fixedImage = IRM;
    m_movingImage = US;
}
    
    

template <typename TFixedImage, typename TMovingImage>
void
    LC2MVImageToImageMetric<TFixedImage,TMovingImage>::SetFixed(typename TFixedImage::Pointer IRM)
    {
        m_fixedImage = IRM;
    }

template <typename TFixedImage, typename TMovingImage>
void
    LC2MVImageToImageMetric<TFixedImage,TMovingImage>::SetMoving(typename TMovingImage::Pointer US)
    {
        m_movingImage = US;
    }
    
template <typename TFixedImage, typename TMovingImage>
void
LC2MVImageToImageMetric<TFixedImage,TMovingImage>::ComputeGradImage()
    {
        std::cout<<"compute gradient image of fixed MRI image"<<std::endl;
        
        typename GradientFilterType::Pointer filterG = GradientFilterType::New();
        filterG->SetInput(m_fixedImage);
        try {
            filterG->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while computing gradient image"<<std::endl;
            std::cerr<<e<<std::endl;
            EXIT_FAILURE;
        }
        
        m_grad = filterG->GetOutput();
        
        //write image for test
        
                typename WriterType::Pointer writer1 = WriterType::New();
                string out1 = "/Users/maximegerard/Documents/testgradMV.nii.gz";
                NiftiImageIO::Pointer io = NiftiImageIO::New();
                writer1->SetInput(m_grad);
                writer1->SetImageIO(io);
                writer1->SetFileName(out1);
                try {
                    writer1->Update();
                } catch (itk::ExceptionObject &e) {
                    cerr<<"error while writing image file"<<endl;
                    cerr<<e<<endl;
                    EXIT_FAILURE;
                }
        
                cout<<"done writing gradient image"<<endl;

    }

template <typename TFixedImage, typename TMovingImage>
void
    LC2MVImageToImageMetric<TFixedImage,TMovingImage>::ComputeMask()
    {
        std::cout<<"computing cropping mask"<<std::endl;
        //binarisation
        
        typename BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();
        thresholder->SetInput(m_movingImage);
        thresholder->SetOutsideValue(255);
        thresholder->SetInsideValue(0);
        thresholder->SetLowerThreshold(0);
        thresholder->SetUpperThreshold(1);
        
        try {
            thresholder->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while binarizing US image"<<std::endl;
            std::cerr<<e<<std::endl;
        }
        
        MaskType::Pointer mask1 =thresholder->GetOutput() ;
        
        
        std::cout<<"done writing initial mask image"<<std::endl;
        
        std::cout<<"test closing"<<std::endl;
        //operation morphologiques
        
        
        kernelType::RadiusType radius;
        radius.Fill(3);
        
        kernelType kernel = kernelType::Ball(radius);
        std::cout<<"radius kernel : "<<kernel.GetRadius()<<std::endl;
        std::cout<<"kernel size : "<<kernel.GetSize()<<std::endl;
        
        CloserType::Pointer closer = CloserType::New();
        closer->SetInput(mask1);
        closer->SetKernel(kernel);
        closer->Update();
        
        m_mask = closer->GetOutput();
        
        //writing mask images
        
        BinaryWriterType::Pointer writer3 = BinaryWriterType::New();
        std::string out3 = "/Users/maximegerard/Documents/testmaskMV.nii.gz";
        itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
        writer3->SetInput(m_mask);
        writer3->SetImageIO(io);
        writer3->SetFileName(out3);
        try {
            writer3->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while writing image file"<<std::endl;
            std::cerr<<e<<std::endl;
            EXIT_FAILURE;
        }
        //
        std::cout<<"done writing final mask image"<<std::endl;
        

    }
    
/**
 * Get the match Measure
 */
template <typename TFixedImage, typename TMovingImage>
typename LC2MVImageToImageMetric<TFixedImage, TMovingImage>::MeasureType
LC2MVImageToImageMetric<TFixedImage, TMovingImage>
::GetValue(const TransformParametersType & parameters) const
{
    //setting des images
    
    //resetting the varsum and lc2sum for new computation
    double variancesum2 =0;
    double lc2varsum2 = 0;
    
    
    if( !m_fixedImage )
    {
        itkExceptionMacro(<< "Fixed image has not been assigned");
    }
    
    
    
    //image qui bouge = image US
    
    if(!m_movingImage)
    {
        itkExceptionMacro(<<"Moving image has not been assigned");
    }
    
    //transformation with regard to TransformParameters
    
    
    //m_Transform = EulerTransformType::New();
    m_Transform->SetParameters(parameters);
    std::cout<<"euler tsf parameters : "<<m_Transform->GetParameters()<<std::endl;
    
    typename TMovingImage::SizeType sizeUS = m_movingImage->GetLargestPossibleRegion().GetSize();
    typename TMovingImage::PointType origin = m_movingImage->GetOrigin();
    typename TMovingImage::SpacingType spacing = m_movingImage->GetSpacing();
    typename TMovingImage::PointType center;
    center[0] = origin[0]+spacing[0]*sizeUS[0]/2;
    center[1] = origin[1]+spacing[1]*sizeUS[1]/2;
    center[2] = origin[2]+spacing[2]*sizeUS[2]/2;
    
    
    EulerTransformType::ParametersType eulerFixedParameters(3);
    eulerFixedParameters[0] =center[0];
    eulerFixedParameters[1] =center[1];
    eulerFixedParameters[2] =center[2];
    
    m_Transform->SetFixedParameters(eulerFixedParameters);
    std::cout<<"tsf fixed param : "<<m_Transform->GetFixedParameters()<<std::endl;
    
    
    
    typename ResamplerType::Pointer resamplefilter = ResamplerType::New();
    resamplefilter->SetInput(m_movingImage);
    resamplefilter->SetSize(m_fixedImage->GetLargestPossibleRegion().GetSize());
    resamplefilter->SetOutputSpacing(m_fixedImage->GetSpacing());
    resamplefilter->SetOutputDirection(m_fixedImage->GetDirection());
    resamplefilter->SetOutputOrigin(m_fixedImage->GetOrigin());
    resamplefilter->SetTransform(m_Transform);
    //resamplefilter->SetTransform(transform);
    
    try {
        resamplefilter->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while transforming moving image"<<std::endl;
        std::cerr<<e<<std::endl;
    }
    
    typename TMovingImage::Pointer movedImage = resamplefilter->GetOutput();
    
    //downsampling de l'image US
    
    typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
    shrinkFilter->SetInput(movedImage);
    shrinkFilter->SetShrinkFactor(0,2);
    shrinkFilter->SetShrinkFactor(1,2);
    shrinkFilter->SetShrinkFactor(2,2);
    try {
        shrinkFilter->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while downsampling transormed us"<<std::endl;
        std::cerr<<e<<std::endl;
    }
    
    //downsampled transformed US = the one on which we effectuate the LC2 computation
    typename TMovingImage::Pointer movingImageT = shrinkFilter->GetOutput();
    
    //Transformation de l'image binaire
    //la tsf est la mm que pour l'US
    
    //transformation du mask
    ResamplerBinaryType::Pointer maskResampler = ResamplerBinaryType::New();
    maskResampler->SetInput(m_mask);
    maskResampler->SetOutputDirection(m_fixedImage->GetDirection());
    maskResampler->SetOutputOrigin(m_fixedImage->GetOrigin());
    maskResampler->SetOutputSpacing(m_fixedImage->GetSpacing());
    maskResampler->SetSize(m_fixedImage->GetLargestPossibleRegion().GetSize());
    maskResampler->SetTransform(m_Transform);
    
    try {
        maskResampler->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while translating image"<<std::endl;
        std::cerr<<e<<std::endl;
    }
    
    MaskType::Pointer tsfMask = maskResampler->GetOutput();
    
    //down sampling du mask
    
    BinaryShrinkFilterType::Pointer binaryShrink = BinaryShrinkFilterType::New();
    binaryShrink->SetInput(tsfMask);
    binaryShrink->SetShrinkFactor(0, 2);
    binaryShrink->SetShrinkFactor(1, 2);
    binaryShrink->SetShrinkFactor(2, 2);
    try {
        binaryShrink->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while downsampling US image"<<std::endl;
        std::cerr<<e<<std::endl;
    }
    
    MaskType::Pointer mask_shrunk = binaryShrink->GetOutput();
    
    //Cropping considerant le zone non nulle de l'US
    //pour garder les boundaries
    int indMaxX =0;
    int indMinX = 100000;
    int indMaxY =0;
    int indMinY = 100000;
    int indMaxZ =0;
    int indMinZ = 100000;
    
    //iterateur pour le mask
    
    BinaryImageIteratorType mask_it(mask_shrunk,mask_shrunk->GetLargestPossibleRegion());
    mask_it.GoToBegin();
    while(!mask_it.IsAtEnd())
    {
        if(mask_it.Get()>0)
        {
            MaskType::IndexType ind = mask_it.GetIndex();
            
            //X
            if(ind[0]>indMaxX)
            {
                indMaxX=ind[0];
            }
            
            else if(ind[0]<indMinX)
            {
                indMinX=ind[0];
            }
            
            //Y
            if(ind[1]>indMaxY)
            {
                indMaxY=ind[1];
            }
            
            else if(ind[1]<indMinY)
            {
                indMinY=ind[1];
            }
            
            //Z
            if(ind[2]>indMaxZ)
            {
                indMaxZ=ind[2];
            }
            
            else if(ind[2]<indMinZ)
            {
                indMinZ=ind[2];
            }
        }
        
        ++mask_it;
        
    }
    
    std::cout<<"indices minimum X,Y,Z "<<std::endl;
    std::cout<<"X : "<<indMinX<<" "<<indMaxX<<std::endl;
    std::cout<<"Y : "<<indMinY<<" "<<indMaxY<<std::endl;
    std::cout<<"Z : "<<indMinZ<<" "<<indMaxZ<<std::endl;
    
    typename TMovingImage::IndexType startCropped;
    startCropped[0] = indMinX ;
    startCropped[1] = indMinY;
    startCropped[2] = indMinZ;
    
    
    
    typename TMovingImage::IndexType endCropped;
    endCropped[0] = indMaxX;
    endCropped[1] = indMaxY;
    endCropped[2] = indMaxZ;
    
    
    
    typename TMovingImage::SizeType sizeCropped;
    sizeCropped[0] = endCropped[0]-startCropped[0]+1;
    sizeCropped[1] = endCropped[1]-startCropped[1]+1;
    sizeCropped[2] = endCropped[2]-startCropped[2]+1;
    
    typename TMovingImage::RegionType regionCropped;
    regionCropped.SetIndex(startCropped);
    regionCropped.SetSize(sizeCropped);
    
    //Cropping de l'image originale
    
    typename ExtractorType::Pointer CroppingFilter = ExtractorType::New();
    CroppingFilter->SetExtractionRegion(regionCropped);
    CroppingFilter->SetInput(movingImageT);
    CroppingFilter->SetDirectionCollapseToIdentity();
    CroppingFilter->Update();
    typename TMovingImage::Pointer moving_Cropped = CroppingFilter->GetOutput();
    std::cout<<"verification image size : "<<moving_Cropped->GetLargestPossibleRegion().GetSize()<<std::endl;
    
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
    
    
    
    //on ne considere que les pixels non nuls de l'image US
    
    // On peut parcourir l'image avec le voisinnage mais ne faire les actions qui si le pixel eq dans le mask est a 1
    
    //////////////////////////////////////////
    // REGION ACCESSIBLE POUR NEIGHBORHOOD //
    ////////////////////////////////////////
    std::cout<<"defining accessible region of image"<<std::endl;
    //on itere sur toute l'image accessible -> celle pour laquelle le neighboorhood it ne va pas sortir de l'image
    typename TMovingImage::RegionType accessibleImagePart;
    typename TMovingImage::RegionType::IndexType startIndex = moving_Cropped->GetLargestPossibleRegion().GetIndex();
    startIndex[0] =startIndex[0]+3;
    startIndex[1] =startIndex[1]+3;
    startIndex[2] =startIndex[2]+3;
    typename TMovingImage::RegionType::SizeType sizeAccessible;
    typename TMovingImage::SizeType sizeIm = moving_Cropped->GetLargestPossibleRegion().GetSize();
    sizeAccessible[0] = sizeIm[0]-6;
    sizeAccessible[1] = sizeIm[1]-6;
    sizeAccessible[2] = sizeIm[2]-6;
    
    accessibleImagePart.SetIndex(startIndex);
    accessibleImagePart.SetSize(sizeAccessible);
    std::cout<<" start index region accessible : "<<startIndex<<std::endl;
    std::cout<<"taille region accessible : "<<accessibleImagePart.GetSize()<<std::endl;
    
    
    //TEST IMAGE ITERATOR RATHER THAN NITERATOR
    //CONST IT CAUSE WE DONT MODIFY IMAGE INTENSITIES
    
    
    //////////////////////////
    // ITERATIONS OVER IMAGE /
    /////////////////////////
    //to iterate over the accessible part of fixedImage
    ImageConstIteratorType US_it(movingImageT,accessibleImagePart);
    US_it.GoToBegin();
    
    
    std::cout<<"calcul de la lc2"<<std::endl;
    
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
    
    std::vector<double> LC2_i;
    LC2_i.clear();
    while(!US_it.IsAtEnd())
    {
        //        //cout<<"US index under consideration : "<<US_it.GetIndex()<<endl;
        //        //on ne considere le voisinage que si le centre appartient à la region blanche du mask et s'il est a l'int de l'im IRM
        //        typename TMovingImage::PointType p;
        //        //get the equivalent in physical point to evaluate whether it's within the MRI image
        //        fixedImage->TransformIndexToPhysicalPoint(US_it.GetIndex(),p);
        //        //cout<<"Physical space coordinates of center of neighbourhood : "<<p<<endl;
        //        typename TMovingImage::IndexType i;
        //        //cout<<"verification mask -US : "<<int(m_mask->GetPixel(US_it.GetIndex()))<<endl;
        
        //we consider the neighbourhood only if it's center is in the actual US data and not outside of the MRI volume
        if(int(mask_shrunk->GetPixel(US_it.GetIndex()))==255)
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
            typename TMovingImage::RegionType neighbourhood;
            typename TMovingImage::RegionType::IndexType start;
            //le debut de la region = le premier indice du masque
            start = US_it.GetIndex();
            //le vrai debut est 3 pixel plus haut, plus a gauche et plus en profondeur
            start[0]= start[0]-3;
            start[1]= start[1]-3;
            start[2]= start[2]-3;
            
            //7-by-7 cube
            typename TMovingImage::RegionType::SizeType sizeN;
            sizeN[0] = 7;
            sizeN[1] = 7;
            sizeN[2] = 7;
            
            neighbourhood.SetIndex(start);
            neighbourhood.SetSize(sizeN);
            
            ImageConstIteratorType it(movingImageT,neighbourhood);
            
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
                typename TMovingImage::IndexType indexUS = it.GetIndex();
                typename TMovingImage::PointType pt;
                
                movingImageT->TransformIndexToPhysicalPoint(indexUS, pt);
                //pt now contains the position in physica space of the considered voxel
                
                typename TFixedImage::IndexType indexIRM;
                
                //si le point est dans l'image IRM
                if(m_fixedImage->TransformPhysicalPointToIndex(pt, indexIRM))
                {
                    M(indice,0) = m_fixedImage->GetPixel(indexIRM);
                    //test vesselness rather than gradient
                    //M(indice,1) = m_vesselness->GetPixel(indexIRM);
                    M(indice,1) = m_grad->GetPixel(indexIRM);
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
                std::cerr<<"Matrix det is null"<<std::endl;
                std::cerr<<e<<std::endl;
                // cerr<<"spatial position of center of neighbourhood"<<p<<endl;
                std::cerr<<"matric M"<<M<<std::endl;
                std::cerr<<"metrice MTM"<<MTM<<std::endl;
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
                //lc2 = 1 - (sum/(m*variance));
                lc2 = sum/(m*variance); //LM is minimizing only !
                
            }
            else
            {
                std::cout<<"variance on patch is null"<<std::endl;
                //cout<<"neighbourhood voxel center : "<<p<<endl;
                std::cout<<U<<std::endl;
                lc2 = 0;
                
            }
            
            //ajout dans vector LC2_i
            LC2_i.push_back(lc2);
            
            
            //cout<<"lc2 locale : "<<lc2<<endl;
            
            //ajout pondere par variance pour calcul de la LC2 totale
            lc2varsum2 = lc2varsum2 + (lc2*variance);
            
            
        }
        
        ++US_it;
        
    }
    
    std::cout<<"done parcours image US"<<std::endl;
    tend = std::time(NULL);
    texec = std::difftime(tend,tbegin);
    
    std::cout<<"temps de parcours en s : "<<texec<<std::endl;
    
    //lc2 finale = moyenne ponderee
    
    
    double lc2final = lc2varsum2/variancesum2;
    std::cout<<"lc2 globale : "<<lc2final<<std::endl;
    
    //placement des element du vecteur dans l'array et return l'array
    //itk::Array<double> measure;
    
    MeasureType measure;
    
    measure.SetSize(LC2_i.size());
    for(int i =0; i<LC2_i.size();i++)
    {
        measure.SetElement(i,LC2_i[i]);
    }
    
    //SetNumberOfValues(LC2_i.size());
    
    
    
    return measure;

}

/**
 * Get the Derivative Measure
 */
template <typename TFixedImage, typename TMovingImage>
void
LC2MVImageToImageMetric<TFixedImage, TMovingImage>
::GetDerivative(const TransformParametersType & parameters,
                DerivativeType & derivative) const
{
    
}

template<typename TFixedImage, typename TMovingImage>
    void
    LC2MVImageToImageMetric<TFixedImage,TMovingImage>::SetTransformParameters(const ParametersType &parameters)const
    {
        m_Transform->SetParameters(parameters);
    }

template<typename TFixedImage, typename TMovingImage>
unsigned int
    LC2MVImageToImageMetric<TFixedImage,TMovingImage>::GetNumberOfValues() const
{
    return m_numberofValues;
        
}
    
template<typename TFixedImage, typename TMovingImage>
    void
    LC2MVImageToImageMetric<TFixedImage,TMovingImage>::SetTransform(const TransformPointer &transform) const
    {
        m_Transform = transform;
    }

template<typename TFixedImage, typename TMovingImage>
    void
    LC2MVImageToImageMetric<TFixedImage,TMovingImage>::SetNumberOfValues(const int number) const
    {
        m_numberofValues = number;
    }

template <typename TFixedImage, typename TMovingImage>
void
LC2MVImageToImageMetric<TFixedImage, TMovingImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
 
}

} // end namespace itk

#endif
