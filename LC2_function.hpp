//
//  LC2_function.hpp
//  LC2_M
//
//  Created by Maxime Gérard on 10/02/16.
//  Copyright © 2016 Maxime Gérard. All rights reserved.
//

#ifndef LC2_function_hpp
#define LC2_function_hpp

#include <stdio.h>
#include <iostream>
#include <string>
#include <time.h>

//ITK CLASSES
#include "itkImage.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "gradient.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkNiftiImageIO.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkExtractImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include <time.h>
#include <vector>
#include "itkResampleImageFilter.h"
#include "itkEuler3DTransform.h"
#include "itkShrinkImageFilter.h"


//DLIB ELEMENTS
#include <dlib/matrix.h>

using namespace dlib;
using namespace std;

//images

typedef itk::Image<double,3> ImageType;
typedef itk::Image<unsigned char,3> MaskType;

//matrices and vectors
typedef itk::Matrix<double, 343,3> MType;
typedef itk::Matrix<double,3,343> MTType;
typedef itk::Matrix<double,3,3> M3Type;
typedef itk::Vector<double, 343> VType;
typedef itk::Vector<double,3> PType;

//gradient
typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType> GradientFilterType;
typedef itk::ImageFileWriter<ImageType> WriterType;

//mask for MRI
typedef itk::Image<unsigned char, 3> MaskType;
typedef itk::ImageFileWriter<MaskType> BinaryWriterType;
typedef itk::BinaryThresholdImageFilter<ImageType, MaskType> BinaryThresholdFilterType;
typedef itk::FlatStructuringElement<3> kernelType;
typedef itk::BinaryMorphologicalClosingImageFilter<MaskType, MaskType, kernelType> CloserType;


//cropping
typedef itk::ExtractImageFilter<ImageType, ImageType> ExtractorType;
typedef itk::ExtractImageFilter<MaskType, MaskType> BinaryExtractorType;


//image iterator
typedef itk::ImageRegionConstIterator<ImageType> ImageConstIteratorType;
typedef itk::ImageRegionIterator<ImageType> ImageIteratorType;
typedef itk::ImageRegionIterator<MaskType> BinaryImageIteratorType;

//for image tsf
typedef itk::Euler3DTransform<double> EulerTransformType;
typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
typedef itk::ResampleImageFilter<MaskType, MaskType> ResamplerBinaryType;
typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkFilterType;
typedef itk::ShrinkImageFilter<MaskType, MaskType> BinaryShrinkFilterType;


class LC2_function
{
public:
    LC2_function(ImageType::Pointer im_Fixed, ImageType::Pointer im_Moving);
    
    double operator()(const dlib::matrix<double>& params)
    const
    {
        //setting des images
        //std::cout<<"iteration : "<<track<<std::endl;
        
        //resetting the varsum and lc2sum for new computation
        double variancesum2 =0;
        double lc2varsum2 = 0;
        
        //defining images
        
        if( !m_FixedImage )
        {
            std::cout<< "Fixed image has not been assigned"<<std::endl;
        }
        
        
        
        //image qui bouge = image US
        
        
        if(!m_MovingImage)
        {
            std::cout<<"Moving image has not been assigned"<<std::endl;
        }
        
        //transformation with regard to TransformParameters
        
        //cout<<"verif param : "<<params<<endl;
        
        EulerTransformType::Pointer transform = EulerTransformType::New();
        EulerTransformType::ParametersType parameters(6);
        //mise a l'echelle des parametres
        parameters[0] = params(0)*5/9;
        parameters[1] = params(1)*5/9;
        parameters[2] = params(2)*5/9;
        parameters[3] = params(3)*100/9;
        parameters[4] = params(4)*100/9;
        parameters[5] = params(5)*100/9;
        transform->SetParameters(parameters);
        std::cout<<"euler tsf parameters : "<<transform->GetParameters()<<std::endl;
        
        typename ImageType::SizeType sizeUS = m_MovingImage->GetLargestPossibleRegion().GetSize();
        typename ImageType::PointType origin = m_MovingImage->GetOrigin();
        typename ImageType::SpacingType spacing = m_MovingImage->GetSpacing();
        typename ImageType::PointType center;
        center[0] = origin[0]+spacing[0]*sizeUS[0]/2;
        center[1] = origin[1]+spacing[1]*sizeUS[1]/2;
        center[2] = origin[2]+spacing[2]*sizeUS[2]/2;
        
        
        EulerTransformType::ParametersType eulerFixedParameters(3);
        eulerFixedParameters[0] =center[0];
        eulerFixedParameters[1] =center[1];
        eulerFixedParameters[2] =center[2];
        
        transform->SetFixedParameters(eulerFixedParameters);
        //std::cout<<"tsf fixed param : "<<transform->GetFixedParameters()<<std::endl;
        
        
        
        typename ResamplerType::Pointer resamplefilter = ResamplerType::New();
        resamplefilter->SetInput(m_MovingImage);
        resamplefilter->SetSize(m_FixedImage->GetLargestPossibleRegion().GetSize());
        resamplefilter->SetOutputSpacing(m_FixedImage->GetSpacing());
        resamplefilter->SetOutputDirection(m_FixedImage->GetDirection());
        resamplefilter->SetOutputOrigin(m_FixedImage->GetOrigin());
        resamplefilter->SetTransform(transform);
        //resamplefilter->SetTransform(transform);
        
        try {
            resamplefilter->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while transforming moving image"<<std::endl;
            std::cerr<<e<<std::endl;
            return EXIT_FAILURE;
        }
        
        typename ImageType::Pointer movedImage = resamplefilter->GetOutput();
        
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
            return EXIT_FAILURE;
        }
        
        //downsampled transformed US = the one on which we effectuate the LC2 computation
        typename ImageType::Pointer movingImageT = shrinkFilter->GetOutput();
        
        //Transformation de l'image binaire
        //la tsf est la mm que pour l'US
        
        //transformation du mask
        ResamplerBinaryType::Pointer maskResampler = ResamplerBinaryType::New();
        maskResampler->SetInput(m_mask);
        maskResampler->SetOutputDirection(m_FixedImage->GetDirection());
        maskResampler->SetOutputOrigin(m_FixedImage->GetOrigin());
        maskResampler->SetOutputSpacing(m_FixedImage->GetSpacing());
        maskResampler->SetSize(m_FixedImage->GetLargestPossibleRegion().GetSize());
        maskResampler->SetTransform(transform);
        
        try {
            maskResampler->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while translating image"<<std::endl;
            std::cerr<<e<<std::endl;
            return EXIT_FAILURE;
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
            return EXIT_FAILURE;
        }
        
        MaskType::Pointer mask_shrunk = binaryShrink->GetOutput();
        
        //Cropping considerant le zone non nulle de l'US
        //pour garder les boundaries
        long indMaxX =0;
        long indMinX = 100000;
        long indMaxY =0;
        long indMinY = 100000;
        long indMaxZ =0;
        long indMinZ = 100000;
        
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
        
        typename ImageType::IndexType startCropped;
        startCropped[0] = indMinX ;
        startCropped[1] = indMinY;
        startCropped[2] = indMinZ;
        
        
        
        typename ImageType::IndexType endCropped;
        endCropped[0] = indMaxX;
        endCropped[1] = indMaxY;
        endCropped[2] = indMaxZ;
        
        
        
        typename ImageType::SizeType sizeCropped;
        sizeCropped[0] = endCropped[0]-startCropped[0]+1;
        sizeCropped[1] = endCropped[1]-startCropped[1]+1;
        sizeCropped[2] = endCropped[2]-startCropped[2]+1;
        
        typename ImageType::RegionType regionCropped;
        regionCropped.SetIndex(startCropped);
        regionCropped.SetSize(sizeCropped);
        
        //Cropping de l'image originale
        
        typename ExtractorType::Pointer CroppingFilter = ExtractorType::New();
        CroppingFilter->SetExtractionRegion(regionCropped);
        CroppingFilter->SetInput(movingImageT);
        CroppingFilter->SetDirectionCollapseToIdentity();
        CroppingFilter->Update();
        typename ImageType::Pointer moving_Cropped = CroppingFilter->GetOutput();
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
        typename ImageType::RegionType accessibleImagePart;
        typename ImageType::RegionType::IndexType startIndex = moving_Cropped->GetLargestPossibleRegion().GetIndex();
        startIndex[0] =startIndex[0]+3;
        startIndex[1] =startIndex[1]+3;
        startIndex[2] =startIndex[2]+3;
        typename ImageType::RegionType::SizeType sizeAccessible;
        typename ImageType::SizeType sizeIm = moving_Cropped->GetLargestPossibleRegion().GetSize();
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
        //to iterate over the accessible part of m_FixedImage
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
        
        
        while(!US_it.IsAtEnd())
        {
            //        //cout<<"US index under consideration : "<<US_it.GetIndex()<<endl;
            //        //on ne considere le voisinage que si le centre appartient à la region blanche du mask et s'il est a l'int de l'im IRM
            //        typename ImageType::PointType p;
            //        //get the equivalent in physical point to evaluate whether it's within the MRI image
            //        m_FixedImage->TransformIndexToPhysicalPoint(US_it.GetIndex(),p);
            //        //cout<<"Physical space coordinates of center of neighbourhood : "<<p<<endl;
            //        typename ImageType::IndexType i;
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
                typename ImageType::RegionType neighbourhood;
                typename ImageType::RegionType::IndexType start;
                //le debut de la region = le premier indice du masque
                start = US_it.GetIndex();
                //le vrai debut est 3 pixel plus haut, plus a gauche et plus en profondeur
                start[0]= start[0]-3;
                start[1]= start[1]-3;
                start[2]= start[2]-3;
                
                //7-by-7 cube
                typename ImageType::RegionType::SizeType sizeN;
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
                    typename ImageType::IndexType indexUS = it.GetIndex();
                    typename ImageType::PointType pt;
                    
                    movingImageT->TransformIndexToPhysicalPoint(indexUS, pt);
                    //pt now contains the position in physica space of the considered voxel
                    
                    typename ImageType::IndexType indexIRM;
                    
                    //si le point est dans l'image IRM
                    if(m_FixedImage->TransformPhysicalPointToIndex(pt, indexIRM))
                    {
                        M(indice,0) = m_FixedImage->GetPixel(indexIRM);
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
                    lc2 = 1 - (sum/(m*variance));
                    //lc2 = sum/(m*variance); //for minimisation test
                    
                }
                else
                {
                    std::cout<<"variance on patch is null"<<std::endl;
                    //cout<<"neighbourhood voxel center : "<<p<<endl;
                    std::cout<<U<<std::endl;
                    lc2 = 0;
                    
                }
                
                
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
        std::cout << std::endl << std::endl;
        
        
        
        
        
        double measure = lc2final;
        
        return measure;
    }
    
    void computeMask();
    void computeGradient();
    int track;
    
private:
    ImageType::Pointer m_MovingImage;
    ImageType::Pointer m_FixedImage;
    MaskType::Pointer m_mask;
    ImageType::Pointer m_grad;
    
};

#endif /* LC2_function_hpp */
