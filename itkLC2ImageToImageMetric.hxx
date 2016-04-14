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
    this->m_MovingImage = US;
    this->m_FixedImage = IRM;
    
}
    

template<typename TFixedImage, typename TMovingImage>
    void LC2ImageToImageMetric<TFixedImage,TMovingImage>::SetFixed(typename TFixedImage::Pointer IRM)
    {
        this->m_FixedImage = IRM;
       
    }
    
template<typename TFixedImage, typename TMovingImage>
    void LC2ImageToImageMetric<TFixedImage,TMovingImage>::SetMoving(typename TMovingImage::Pointer US)
    {
        this->m_MovingImage =US;
    }
    
    
/********************
 * IMAGES TSF
 ********************/
    
    template<typename TFixedImage, typename TMovingImage>
    typename TMovingImage::Pointer
    LC2ImageToImageMetric<TFixedImage,TMovingImage>::TransformImage(const TransformParametersType &parameters, int ind) const
    {
        ImageType::Pointer imageTransformed;
        
        //1 represents euler tsf
        if(ind==1)
        {
            EulerTransformType::Pointer transform = EulerTransformType::New();
            transform->SetParameters(parameters);
            //std::cout<<"euler tsf parameters : "<<transform->GetParameters()<<std::endl;
            
            typename ImageType::SizeType sizeUS = this->m_MovingImage->GetLargestPossibleRegion().GetSize();
            typename ImageType::PointType origin = this->m_MovingImage->GetOrigin();
            typename ImageType::SpacingType spacing = this->m_MovingImage->GetSpacing();
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
            resamplefilter->SetInput(this->m_MovingImage);
            resamplefilter->SetSize(this->m_FixedImage->GetLargestPossibleRegion().GetSize());
            resamplefilter->SetOutputSpacing(this->m_FixedImage->GetSpacing());
            resamplefilter->SetOutputDirection(this->m_FixedImage->GetDirection());
            resamplefilter->SetOutputOrigin(this->m_FixedImage->GetOrigin());
            resamplefilter->SetTransform(transform);
            //resamplefilter->SetTransform(transform);
            
            try {
                resamplefilter->Update();
            } catch (itk::ExceptionObject &e) {
                std::cerr<<"error while transforming moving image"<<std::endl;
                std::cerr<<e<<std::endl;
            }
            
            imageTransformed= resamplefilter->GetOutput();
            
            
        }
        
        return imageTransformed;

    }
    
    
    
    template<typename TFixedImage, typename TMovingImage>
    typename MaskType::Pointer
    LC2ImageToImageMetric<TFixedImage,TMovingImage>::TransformMask(const TransformParametersType &parameters, int ind) const
    {
        MaskType::Pointer maskTsf;
        
        //1 represents euler tsf
        if(ind==1)
        {
            EulerTransformType::Pointer transform = EulerTransformType::New();
            transform->SetParameters(parameters);
            //std::cout<<"euler tsf parameters : "<<transform->GetParameters()<<std::endl;
            
            typename ImageType::SizeType sizeUS = this->m_MovingImage->GetLargestPossibleRegion().GetSize();
            typename ImageType::PointType origin = this->m_MovingImage->GetOrigin();
            typename ImageType::SpacingType spacing = this->m_MovingImage->GetSpacing();
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
            
            ResamplerBinaryType::Pointer maskResampler = ResamplerBinaryType::New();
            maskResampler->SetInput(m_mask);
            maskResampler->SetOutputDirection(this->m_FixedImage->GetDirection());
            maskResampler->SetOutputOrigin(this->m_FixedImage->GetOrigin());
            maskResampler->SetOutputSpacing(this->m_FixedImage->GetSpacing());
            maskResampler->SetSize(this->m_FixedImage->GetLargestPossibleRegion().GetSize());
            maskResampler->SetTransform(transform);
            
            try {
                maskResampler->Update();
            } catch (itk::ExceptionObject &e) {
                cerr<<"error while transforming mesh"<<endl;
                cerr<<e<<endl;
            }
            
            maskTsf = maskResampler->GetOutput();
            
        }
        
        return maskTsf;
            
        }
    
    template<typename TFixedImage, typename TMovingImage>
    void
    LC2ImageToImageMetric<TFixedImage,TMovingImage>::setLiverMask(MaskType::Pointer Mask)
    {
        m_useLiverMask=true;
        m_LiverMask = Mask;
    }
    

/************************
 * GRADIENT COMPUTATION
 ************************/
    
    template<typename TFixedImage, typename TMovingImage>
    void
    LC2ImageToImageMetric<TFixedImage,TMovingImage>::ComputeGradImage()
    {
        std::cout<<"compute gradient image of fixed MRI image"<<std::endl;
        
        typename GradientFilterType::Pointer filterG = GradientFilterType::New();
        filterG->SetInput(this->m_FixedImage);
        try {
            filterG->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while computing gradient image"<<std::endl;
            std::cerr<<e<<std::endl;
            EXIT_FAILURE;
        }
        
        m_grad = filterG->GetOutput();
        
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

        
    }
    


/***********************
 * US MASK FOR CROPPING
 ***********************/
    template <typename TFixedImage, typename TMovingImage>
    void
    LC2ImageToImageMetric<TFixedImage, TMovingImage>::ComputeMask()
    {
        std::cout<<"computing cropping mask"<<std::endl;
        //binarisation
        
        typename BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();
        thresholder->SetInput(this->m_MovingImage);
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
        
//        BinaryWriterType::Pointer writer3 = BinaryWriterType::New();
//        std::string out3 = "/Users/maximegerard/Documents/testmask.nii.gz";
//        itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//        writer3->SetInput(m_mask);
//        writer3->SetImageIO(io);
//        writer3->SetFileName(out3);
//        try {
//            writer3->Update();
//        } catch (itk::ExceptionObject &e) {
//            std::cerr<<"error while writing image file"<<std::endl;
//            std::cerr<<e<<std::endl;
//            EXIT_FAILURE;
//        }
////
//        std::cout<<"done writing final mask image"<<std::endl;

        
    }
    
    
/** test the metric function
 */
    
template <typename TFixedImage, typename TMovingImage>
double
    LC2ImageToImageMetric<TFixedImage, TMovingImage>::test(const TransformParametersType &parameters)
    {
        //resetting the varsum and lc2sum for new computation
        variancesum =0;
        lc2varsum = 0;
        
        //defining images
        //image fixe = image IRM
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
        
        std::cout<<"compute gradient of MRI image"<<std::endl;
        
        typename GradientFilterType::Pointer filterG = GradientFilterType::New();
        filterG->SetInput(fixedImage);
        try {
            filterG->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while computing gradient image"<<std::endl;
            std::cerr<<e<<std::endl;
            EXIT_FAILURE;
        }
        
        typename TMovingImage::Pointer image_grad = filterG->GetOutput();
        
//        //write image for test
//        
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

        std::cout<<"done writing gradient image"<<std::endl;
        
//        ///////////////////////////////////////
//        //CROPPING CONSIDERING MRI POSITION //
//        //////////////////////////////////////
//        
//        
//        //on s'interesse au voxel au milieu de chacun des face du volume MRI
//        
//        //get the limits of images in terms of indexes
//        typename TMovingImage::SizeType sizeIRM = fixedImage->GetLargestPossibleRegion().GetSize();
//        cout<<"size IRM : "<<sizeIRM<<endl;
//        typename TFixedImage::SizeType sizeUS = movingImage->GetLargestPossibleRegion().GetSize();
//        cout<<"size US : "<<sizeUS<<endl;
//        
//        typename TFixedImage::IndexType ind1;
//        ind1[0] = 0;
//        ind1[1] = sizeIRM[1]/2;
//        ind1[2] = sizeIRM[2]/2;
//        //cout<<"index IRM p1 : "<<ind1<<endl;
//        typename TFixedImage::PointType p1;
//        fixedImage->TransformIndexToPhysicalPoint(ind1,p1);
//        //cout<<"position spatiale p1 : "<<p1<<endl;;
//        //eq US
//        typename TMovingImage::IndexType indUS1;
//        int indUS_infX;
//        if(movingImage->TransformPhysicalPointToIndex(p1,indUS1))
//        {
//            //cout<<"le point 1 est dans le volume US !"<<endl;
//            //cout<<"indice US : "<<indUS1<<endl;
//            indUS_infX = indUS1[0];
//            
//        }
//        else
//            indUS_infX = 0;
//        
//        typename TFixedImage::IndexType ind2;
//        ind2[0] = sizeIRM[0]-1;
//        ind2[1] = sizeIRM[1]/2;
//        ind2[2] = sizeIRM[2]/2;
//        //cout<<"indice irm p2 : "<<ind2<<endl;
//        typename TFixedImage::PointType p2;
//        fixedImage->TransformIndexToPhysicalPoint(ind2,p2);
//        //cout<<"position spatiale p2 : "<<p2<<endl;
//        //eq US
//        typename TMovingImage::IndexType indUS2;
//        int indUS_supX;
//        if(movingImage->TransformPhysicalPointToIndex(p2,indUS2))
//        {
//            //cout<<"le point 2 est dans le volume US !"<<endl;
//            //cout<<"indice US : "<<indUS2<<endl;
//            indUS_supX = indUS2[0];
//            
//        }
//        
//        else
//            indUS_supX = sizeUS[0]-1;
//
//        
//        
//        typename TFixedImage::IndexType ind3;
//        ind3[0] = sizeIRM[0]/2;
//        ind3[1] = 0;
//        ind3[2] = sizeIRM[2]/2;
//        //cout<<"indice irm p3 : "<<ind3<<endl;
//        typename TFixedImage::PointType p3;
//        fixedImage->TransformIndexToPhysicalPoint(ind3,p3);
//        //cout<<"position spatiale p3 : "<<p3<<endl;
//        //eq US
//        typename TMovingImage::IndexType indUS3;
//        int indUS_infY;
//        if(movingImage->TransformPhysicalPointToIndex(p3,indUS3))
//        {
//            //cout<<"le point 3 est dans le volume US !"<<endl;
//            //cout<<"indice us : "<<indUS3<<endl;
//            indUS_infY = indUS3[2];
//            
//        }
//        else
//            indUS_infY = 0;
//        
//        typename TFixedImage::IndexType ind4;
//        ind4[0] = sizeIRM[0]/2;
//        ind4[1] = sizeIRM[1]-1;
//        ind4[2] = sizeIRM[2]/2;
//        //cout<<"indice irm p4 : "<<ind4<<endl;
//        typename TFixedImage::PointType p4;
//        fixedImage->TransformIndexToPhysicalPoint(ind4,p4);
//        //cout<<"position spatiale p4 : "<<p4<<endl;
//        //eq US
//        typename TMovingImage::IndexType indUS4;
//        int indUS_supY;
//        if(movingImage->TransformPhysicalPointToIndex(p4,indUS4))
//        {
//            //cout<<"le point 4 est dans le volume US !"<<endl;
//            //cout<<"indice us : "<<indUS4<<endl;
//            indUS_supY = indUS4[2];
//            
//        }
//        else
//            indUS_supY = sizeUS[2]-1;
//        
//        typename TFixedImage::IndexType ind5;
//        ind5[0] = sizeIRM[0]/2;
//        ind5[1] = sizeIRM[1]/2;
//        ind5[2] = 0;
//        //cout<<"indice irm p5 : "<<ind5<<endl;
//        typename TFixedImage::PointType p5;
//        fixedImage->TransformIndexToPhysicalPoint(ind5,p5);
//        //cout<<"position spatiale p5 : "<<p5<<endl;
//        //eq US
//        typename TMovingImage::IndexType indUS5;
//        int indUS_infZ;
//        if(movingImage->TransformPhysicalPointToIndex(p5,indUS5))
//        {
//            //cout<<"le point 5 est dans le volume US !"<<endl;
//            //cout<<"indice us : "<<indUS5<<endl;
//            indUS_infZ = indUS5[1];
//            
//        }
//        else
//            indUS_infZ = 0;
//        
//        typename TFixedImage::IndexType ind6;
//        ind6[0] = sizeIRM[0]/2;
//        ind6[1] = sizeIRM[1]/2;
//        ind6[2] = sizeIRM[2]-1;
//        //cout<<"indice irm p6 : "<<ind6<<endl;
//        typename TFixedImage::PointType p6;
//        fixedImage->TransformIndexToPhysicalPoint(ind6,p6);
//        //cout<<"position spatiale p6 : "<<p6<<endl;
//        //eq US
//        typename TMovingImage::IndexType indUS6;
//        int indUS_supZ;
//        if(movingImage->TransformPhysicalPointToIndex(p6,indUS6))
//        {
//            //cout<<"le point 6 est dans le volume US !"<<endl;
//            indUS_supZ = indUS6[1];
//            //cout<<"indice us : "<<indUS6<<endl;
//            
//        }
//        else
//            indUS_supZ = sizeUS[1]-1;
//        
////        cout<<"test verification"<<endl;
////        cout<<"indice inf X "<<indUS_infX<<endl;
////        cout<<"indice sup X "<<indUS_supX<<endl;
////        cout<<"indice inf Y "<<indUS_infY<<endl;
////        cout<<"indice sup Y "<<indUS_supY<<endl;
////        cout<<"indice inf Z "<<indUS_infZ<<endl;
////        cout<<"indice sup Z "<<indUS_supZ<<endl;
//
//        
//        typename TMovingImage::IndexType startCropped;
//        startCropped[0] = indUS_infX ;
//        startCropped[1] = indUS_infZ;
//        startCropped[2] = indUS_infY ;
//        
//
//        
//        typename TMovingImage::IndexType endCropped;
//        endCropped[0] = indUS_supX;
//        endCropped[1] = indUS_supZ;
//        endCropped[2] = indUS_supY;
//        
//        
//        
//        typename TMovingImage::SizeType sizeCropped;
//        sizeCropped[0] = endCropped[0]-startCropped[0]+1;
//        sizeCropped[1] = endCropped[1]-startCropped[1]+1;
//        sizeCropped[2] = endCropped[2]-startCropped[2]+1;
//        
//        typename TMovingImage::RegionType regionCropped;
//        regionCropped.SetIndex(startCropped);
//        regionCropped.SetSize(sizeCropped);
//        
//        //Cropping de l'image originale
//        
//        typename ExtractorType::Pointer CroppingFilter = ExtractorType::New();
//        CroppingFilter->SetExtractionRegion(regionCropped);
//        CroppingFilter->SetInput(movingImage);
//        CroppingFilter->SetDirectionCollapseToIdentity();
//        CroppingFilter->Update();
//        typename TMovingImage::Pointer moving_Cropped = CroppingFilter->GetOutput();
//        cout<<"verification image size : "<<moving_Cropped->GetLargestPossibleRegion().GetSize()<<endl;
//        
//        //writing to verify
//        typename WriterType::Pointer writer4 = WriterType::New();
//        //itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//        string out4 = "/Users/maximegerard/Documents/testCroppingUS.nii.gz";
//        writer4->SetImageIO(io);
//        writer4->SetInput(moving_Cropped);
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
//
//        //on ne considere que les pixels non nuls de l'image US
//        
//        // On peut parcourir l'image avec le voisinnage mais ne faire les actions qui si le pixel eq dans le mask est a 1
        

        
        //Transform Mask
        
        EulerTransformType::Pointer tsf = EulerTransformType::New();
        tsf->SetParameters(parameters);
        
        MaskType::SizeType sizeUS = m_mask->GetLargestPossibleRegion().GetSize();
        MaskType::PointType origin = m_mask->GetOrigin();
        MaskType::SpacingType spacing = m_mask->GetSpacing();
        MaskType::PointType center;
        center[0] = origin[0]+spacing[0]*sizeUS[0]/2;
        center[1] = origin[1]+spacing[1]*sizeUS[1]/2;
        center[2] = origin[2]+spacing[2]*sizeUS[2]/2;
        
        
        EulerTransformType::ParametersType eulerFixedParameters(3);
        eulerFixedParameters[0] =center[0];
        eulerFixedParameters[1] =center[1];
        eulerFixedParameters[2] =center[2];
        
        tsf->SetFixedParameters(eulerFixedParameters);
        
        std::cout<<"euler parameters in test function : "<<std::endl;
        std::cout<<"optimizable param : "<<tsf->GetParameters()<<std::endl;
        std::cout<<" fixed param : "<<tsf->GetFixedParameters()<<std::endl;
        
        //transformation du mask
        ResamplerBinaryType::Pointer maskResampler = ResamplerBinaryType::New();
        maskResampler->SetInput(m_mask);
        maskResampler->SetOutputDirection(fixedImage->GetDirection());
        maskResampler->SetOutputOrigin(fixedImage->GetOrigin());
        maskResampler->SetOutputSpacing(fixedImage->GetSpacing());
        maskResampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
        maskResampler->SetTransform(tsf);
        
        try {
            maskResampler->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while translating image"<<std::endl;
            std::cerr<<e<<std::endl;
            return EXIT_FAILURE;
        }
        
        MaskType::Pointer tsfMask = maskResampler->GetOutput();
        
        
//        BinaryWriterType::Pointer writer7 = BinaryWriterType::New();
//        string out7 ="/Users/maximegerard/Documents/tsfmask.nii.gz";
//        writer7->SetImageIO(io);
//        writer7->SetInput(tsfMask);
//        writer7->SetFileName(out7);
//        try {
//            writer7->Update();
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"error whilte writing tsf mask image"<<endl;
//            cerr<<e<<endl;
//            return EXIT_FAILURE;
//        }
        
        std::cout<<"done writing transformed mask"<<std::endl;
        
        //DOWN SAMPLING LE MASK
        
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
        std::cout<<"taille mask basse res : "<<mask_shrunk->GetLargestPossibleRegion().GetSize()<<std::endl;
        
        
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
        CroppingFilter->SetInput(movingImage);
        CroppingFilter->SetDirectionCollapseToIdentity();
        CroppingFilter->Update();
        typename TMovingImage::Pointer moving_Cropped = CroppingFilter->GetOutput();
        std::cout<<"verification image size : "<<moving_Cropped->GetLargestPossibleRegion().GetSize()<<std::endl;
        
//        //writing to verify
//        typename WriterType::Pointer writer4 = WriterType::New();
//        //itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//        string out4 = "/Users/maximegerard/Documents/testCroppingUS.nii.gz";
//        writer4->SetImageIO(io);
//        writer4->SetInput(moving_Cropped);
//        writer4->SetFileName(out4);
//        try {
//            writer4->Update();
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"error while writing croppped image"<<endl;
//            cerr<<e<<endl;
//        }
        
        std::cout<<"done cropping image"<<std::endl;
        
//        //DOWNSAMPLING US ET MASK
//        ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
//        shrinkFilter->SetInput(moving_Cropped);
//        shrinkFilter->SetShrinkFactor(0, 2);
//        shrinkFilter->SetShrinkFactor(1, 2);
//        shrinkFilter->SetShrinkFactor(2, 2);
//        try {
//            shrinkFilter->Update();
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"error while downsampling US image"<<endl;
//            cerr<<e<<endl;
//            return EXIT_FAILURE;
//        }
//        
//        ImageType::moving_Shrunk
        
        
        
        
       
        //////////////////////////////////////////
        // REGION ACCESSIBLE POUR NEIGHBORHOOD //
        ////////////////////////////////////////
        std::cout<<"defining accessible region of image"<<std::endl;
        //on itere sur toute l'image accessible -> celle pour laquelle le neighboorhood it ne va pas sortir de l'image
        typename TMovingImage::RegionType accessibleImagePart;
        typename TMovingImage::RegionType::IndexType startIndex = moving_Cropped->GetLargestPossibleRegion().GetIndex();//moving_Cropped->GetLargestPossibleRegion().GetIndex();
        startIndex[0] =startIndex[0]+3;
        startIndex[1] =startIndex[1]+3;
        startIndex[2] =startIndex[2]+3;
        typename TMovingImage::RegionType::SizeType sizeAccessible;
        typename TMovingImage::SizeType sizeIm = moving_Cropped->GetLargestPossibleRegion().GetSize();//moving_Cropped
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
        ImageConstIteratorType US_it(movingImage,accessibleImagePart);
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
//        interpolator->SetInputImage(movingImage);
//        
//        typename LinearInterpolatorType::Pointer gradInterpolator = LinearInterpolatorType::New();
//        gradInterpolator->SetInputImage(image_grad);
        
        
        while(!US_it.IsAtEnd())
        {
            //cout<<"US index under consideration : "<<US_it.GetIndex()<<endl;
            //on ne considere le voisinage que si le centre appartient à la region blanche du mask et s'il est a l'int de l'im IRM
//            typename TMovingImage::PointType p;
//            //get the equivalent in physical point to evaluate whether it's within the MRI image
//            movingImage->TransformIndexToPhysicalPoint(US_it.GetIndex(),p);
//            //cout<<"Physical space coordinates of center of neighbourhood : "<<p<<endl;
//            typename TFixedImage::IndexType i;
            //cout<<"verification mask -US : "<<int(m_mask->GetPixel(US_it.GetIndex()))<<endl;
            
            //we consider the neighbourhood only if it's center is in the actual US data and not outside of the MRI volume
            if(int(mask_shrunk->GetPixel(US_it.GetIndex()))==255 )//&& fixedImage->TransformPhysicalPointToIndex(p,i))
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
                
                ImageConstIteratorType it(movingImage,neighbourhood);
                
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
                
                    movingImage->TransformIndexToPhysicalPoint(indexUS, pt);
                    //pt now contains the position in physica space of the considered US voxel
                    
                    typename TFixedImage::IndexType indexIRM;
                
                    //si le point est dans l'image IRM
                    if(fixedImage->TransformPhysicalPointToIndex(pt, indexIRM))
                    {
                        M(indice,0) = fixedImage->GetPixel(indexIRM);
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
                        std::cerr<<"Matrix det is null"<<std::endl;
                        std::cerr<<e<<std::endl;
                        //cerr<<"spatial position of center of neighbourhood"<<p<<endl;
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
                    //lc2 = sum/(m*variance);
                    
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
                lc2varsum = lc2varsum + (lc2*variance);
                
                
                }
            
            
            ++US_it;

        }
        
        std::cout<<"done parcours image US"<<std::endl;
        tend = std::time(NULL);
        texec = std::difftime(tend,tbegin);
        
        std::cout<<"temps de parcours en s : "<<texec<<std::endl;
        
    //lc2 finale = moyenne ponderee
        
   
    m_lc2final = lc2varsum/variancesum;
        std::cout<<"lc2 globale : "<<m_lc2final<<std::endl;

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
    //FixedImageConstPointer fixedImage = this->m_FixedImage;
    
    if( !this->m_FixedImage )
    {
        itkExceptionMacro(<< "Fixed image has not been assigned");
    }
    
    
    
    //image qui bouge = image US
    //MovingImageConstPointer movingImage = this->m_MovingImage;
    
    
    if(!this->m_MovingImage)
    {
        itkExceptionMacro(<<"Moving image has not been assigned");
    }
    
    //transformation with regard to TransformParameters
    
    ImageType::Pointer movedImage = TransformImage(parameters, 1);
    MaskType::Pointer movedMask = TransformMask(parameters,1);
 
    
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
    typename TMovingImage::Pointer movingImageT = shrinkFilter->GetOutput();
    
    //Transformation de l'image binaire
    //la tsf est la mm que pour l'US
    
    
    //down sampling du mask
    
    BinaryShrinkFilterType::Pointer binaryShrink = BinaryShrinkFilterType::New();
    binaryShrink->SetInput(movedMask);
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
    typename LinearInterpolatorFilterType::Pointer interpolator = LinearInterpolatorFilterType::New();
    interpolator->SetInputImage(this->m_FixedImage);
    
    typename LinearInterpolatorFilterType::Pointer gradInterpolator = LinearInterpolatorFilterType::New();
    gradInterpolator->SetInputImage(m_grad);
    
    int counter_det_null =0;
    
    while(!US_it.IsAtEnd())
    {
        //        //cout<<"US index under consideration : "<<US_it.GetIndex()<<endl;
        //        //on ne considere le voisinage que si le centre appartient à la region blanche du mask et s'il est a l'int de l'im IRM accessible
        typename ImageType::PointType p;
        //get the equivalent in physical point to evaluate whether it's within the MRI image
        movingImageT->TransformIndexToPhysicalPoint(US_it.GetIndex(),p);
        //cout<<"Physical space coordinates of center of neighbourhood : "<<p<<endl;
        //typename ImageType::IndexType i;
        //cout<<"verification mask -US : "<<int(m_mask->GetPixel(US_it.GetIndex()))<<endl;
        if(m_useLiverMask)
        {
            //std::cout<<"Use of liver mask to limit ROI"<<std::endl;
            ImageType::IndexType i;
            m_LiverMask->TransformPhysicalPointToIndex(p, i);
            
            //we consider the neighbourhood only if it's center is in the actual US data and inside the liver mask
            if(int(mask_shrunk->GetPixel(US_it.GetIndex()))==255 && int(m_LiverMask->GetPixel(i))==1)// && m_accessibleMRI->TransformPhysicalPointToIndex(p, i))
            {
                //cout<<"neighbourhood in real IRM mask data"<<endl;
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
                    if(this->m_FixedImage->TransformPhysicalPointToIndex(pt, indexIRM))
                    {
                        //                        M(indice,0) = m_FixedImage->GetPixel(indexIRM);
                        //
                        //                        M(indice,1) = m_grad->GetPixel(indexIRM);
                        
                        //test with linear interpolator
                        M(indice,0) = interpolator->Evaluate(pt);
                        M(indice,1) = gradInterpolator->Evaluate(pt);
                        
                    }
                    
                    
                    else//si on essaye de comparer a un element hors de l'image IRM
                    {
                        cout<<"en dehors de l'image IRM !"<<endl;
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
                
                
                M3Type MTM = M3Type::Matrix(M.GetTranspose()*M.GetVnlMatrix());
                
                //inversion
                //c'est ici qu'il faut faire la gestion de det =0 !
                
                try {
                    MTM = MTM.GetInverse();
                    // c = (M^T*M)^-1 * M^T*U
                    
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
                    //calcul de la LCI sur ce patch
                    param = MTM*M.GetTranspose()*U;
                    //STD DEV COMPUTATION
                    //calcul de la variance sur ce patch
                    
                    //WE COMPUTE VARIANCE AND LOCAL LC2 ONLY IF MATRIX IS NOT SINGULAR
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
                    
                    
                } catch (itk::ExceptionObject &e) {
                    std::cerr<<"Matrix det is null"<<std::endl;
                    std::cerr<<e<<std::endl;
                    counter_det_null++;
                    // cerr<<"spatial position of center of neighbourhood"<<p<<endl;
                    //std::cerr<<"matric M"<<M<<std::endl;
                    //std::cerr<<"metrice MTM"<<MTM<<std::endl;
                    //ici c'est pcq ce sont des patchs unifomrme qu'il y a erreur ?
                    //param[0] = 1 ;
                    //param[1] = 0;
                    //param[2] = 0;
                    
                    
                }
                
                
                
                
                
            }
        }
        
        else
        {
            //we consider the neighbourhood only if it's center is in the actual US data and not outside of the MRI volume
            if(int(mask_shrunk->GetPixel(US_it.GetIndex()))==255)// && m_accessibleMRI->TransformPhysicalPointToIndex(p, i))
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
                    if(this->m_FixedImage->TransformPhysicalPointToIndex(pt, indexIRM))
                    {
                        //                        M(indice,0) = m_FixedImage->GetPixel(indexIRM);
                        //
                        //                        M(indice,1) = m_grad->GetPixel(indexIRM);
                        
                        //test with linear interpolator
                        M(indice,0) = interpolator->Evaluate(pt);
                        M(indice,1) = gradInterpolator->Evaluate(pt);
                        
                    }
                    
                    
                    else//si on essaye de comparer a un element hors de l'image IRM
                    {
                        cout<<"en dehors de l'image IRM !"<<endl;
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
                
                
                M3Type MTM = M3Type::Matrix(M.GetTranspose()*M.GetVnlMatrix());
                
                //inversion
                //c'est ici qu'il faut faire la gestion de det =0 !
                
                try {
                    MTM = MTM.GetInverse();
                    // c = (M^T*M)^-1 * M^T*U
                    
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
                    //calcul de la LCI sur ce patch
                    param = MTM*M.GetTranspose()*U;
                    //STD DEV COMPUTATION
                    //calcul de la variance sur ce patch
                    
                    //WE COMPUTE VARIANCE AND LOCAL LC2 ONLY IF MATRIX IS NOT SINGULAR
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
                    
                    
                } catch (itk::ExceptionObject &e) {
                    std::cerr<<"Matrix det is null"<<std::endl;
                    std::cerr<<e<<std::endl;
                    counter_det_null++;
                    // cerr<<"spatial position of center of neighbourhood"<<p<<endl;
                    //std::cerr<<"matric M"<<M<<std::endl;
                    //std::cerr<<"metrice MTM"<<MTM<<std::endl;
                    //ici c'est pcq ce sont des patchs unifomrme qu'il y a erreur ?
                    //param[0] = 1 ;
                    //param[1] = 0;
                    //param[2] = 0;
                    
                    
                }
                
                
                
                
                
            }
            
        }
        
        
        
        ++US_it;
        
    }
    
    std::cout<<"done parcours image US"<<std::endl;
    tend = std::time(NULL);
    texec = std::difftime(tend,tbegin);
    
    std::cout<<"temps de parcours en s : "<<texec<<std::endl;
    
    //lc2 finale = moyenne ponderee
    
    std::cout<<"number of null patches : "<<counter_det_null<<std::endl;
    
    
    double lc2final = lc2varsum2/variancesum2;
    std::cout<<"lc2 globale : "<<lc2final<<std::endl;
    std::cout << std::endl << std::endl;

 
    

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
