//
//  LC2_function.cpp
//  LC2_M
//
//  Created by Maxime Gérard on 10/02/16.
//  Copyright © 2016 Maxime Gérard. All rights reserved.
//

#include "LC2_function.hpp"

LC2_function::LC2_function(ImageType::Pointer im_Fixed, ImageType::Pointer im_Moving,string out)
{
    //std::cout<<"creating measure object"<<std::endl;
    m_outputPath=out;
    m_FixedImage = im_Fixed;
    m_MovingImage = im_Moving;
    computeGradient();
    //limitMRI();
    computeMask();
    
    m_maxTrans =0;
    m_maxRot = 0;
    m_radius = 0;
    
    m_useLiverMask=false;
    
    
}

void LC2_function::setLiverMask(MaskType::Pointer liver)
{
    m_useLiverMask=true;
    m_LiverMask=liver;
    
    //we than transform the fixed MRI with that mask to keep only the liver info
    //in the computation of LC2 we are going to have to exclude all points outside that mask
}

void LC2_function::limitMRI()
{
    std::cout<<"test limitating MRI image to interesting area of the liver"<<endl;
    ImageType::IndexType startPoint;
    startPoint[0] = 119;
    startPoint[1] = 180;
    startPoint[2] = 0;
    
    ImageType::SizeType sizeR;
    sizeR[0] = 80;
    sizeR[1] = 86;
    sizeR[2] = 79;
    
    ImageType::RegionType regionR;
    regionR.SetIndex(startPoint);
    regionR.SetSize(sizeR);
    
    ExtractorType::Pointer extractFilter = ExtractorType::New();
    extractFilter->SetExtractionRegion(regionR);
    extractFilter->SetInput(m_FixedImage);
    extractFilter->SetDirectionCollapseToIdentity();
    try {
        extractFilter->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while cropping MRI image"<<std::endl;
        std::cerr<<e<<std::endl;
    }
    
    m_accessibleMRI = extractFilter->GetOutput();
    
    //enregistrement pour verification
    WriterType::Pointer writer10 = WriterType::New();
    itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
    std::string out10 = m_outputPath+"/testCroppingMRI.nii.gz";
    writer10->SetImageIO(io);
    writer10->SetInput(m_accessibleMRI);
    writer10->SetFileName(out10);
    try {
        writer10->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while writing cropped MRI"<<std::endl;
        std::cerr<<e<<std::endl;
    }
}

void LC2_function::computeGradient()
{
    std::cout<<"compute gradient image of fixed MRI image"<<std::endl;
    
    typename GradientFilterType::Pointer filterG = GradientFilterType::New();
    filterG->SetInput(m_FixedImage);
    try {
        filterG->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while computing gradient image"<<std::endl;
        std::cerr<<e<<std::endl;
        EXIT_FAILURE;
    }
    
    m_grad = filterG->GetOutput();
    
    //write image for test
//    
//            typename WriterType::Pointer writer1 = WriterType::New();
//            string out1 = m_outputPath+"/testgrad.nii.gz";
//            itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//            writer1->SetInput(m_grad);
//            writer1->SetImageIO(io);
//            writer1->SetFileName(out1);
//            try {
//                writer1->Update();
//            } catch (itk::ExceptionObject &e) {
//                cerr<<"error while writing image file"<<endl;
//                cerr<<e<<endl;
//                EXIT_FAILURE;
//            }
//    
//            cout<<"done writing gradient image"<<endl;

    
}

void LC2_function::computeMask()
{
    std::cout<<"computing cropping mask"<<std::endl;
    //binarisation
    
    typename BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();
    thresholder->SetInput(m_MovingImage);
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
    
//    BinaryWriterType::Pointer writer3 = BinaryWriterType::New();
//    std::string out3 = m_outputPath+"/testmask.nii.gz";
//    itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//    writer3->SetInput(m_mask);
//    writer3->SetImageIO(io);
//    writer3->SetFileName(out3);
//    try {
//        writer3->Update();
//    } catch (itk::ExceptionObject &e) {
//        std::cerr<<"error while writing image file"<<std::endl;
//        std::cerr<<e<<std::endl;
//        EXIT_FAILURE;
//    }
//    
//    std::cout<<"done writing final mask image"<<std::endl;

    
}

ImageType::Pointer LC2_function::TransformImage(const dlib::matrix<double> &params, int ind) const
{
    ImageType::Pointer imageTransformed;
    
    //1 represents euler tsf
    if(ind==1)
    {
        EulerTransformType::Pointer transform = EulerTransformType::New();
        EulerTransformType::ParametersType parameters(6);
        //mise a l'echelle des parametres
        parameters[0] = params(0)*(m_maxRot/m_radius);
        parameters[1] = params(1)*(m_maxRot/m_radius);
        parameters[2] = params(2)*(m_maxRot/m_radius);
        parameters[3] = params(3)*(m_maxTrans/m_radius);
        parameters[4] = params(4)*(m_maxTrans/m_radius);
        parameters[5] = params(5)*(m_maxTrans/m_radius);
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
        }
        
       imageTransformed= resamplefilter->GetOutput();

        
    }
    
    return imageTransformed;
}

MaskType::Pointer LC2_function::TransformMask(const dlib::matrix<double> &params, int ind) const
{
    MaskType::Pointer maskTsf;
    
    //1 represents euler tsf
    if(ind==1)
    {
        EulerTransformType::Pointer transform = EulerTransformType::New();
        EulerTransformType::ParametersType parameters(6);
        //mise a l'echelle des parametres
        parameters[0] = params(0)*(m_maxRot/m_radius);
        parameters[1] = params(1)*(m_maxRot/m_radius);
        parameters[2] = params(2)*(m_maxRot/m_radius);
        parameters[3] = params(3)*(m_maxTrans/m_radius);
        parameters[4] = params(4)*(m_maxTrans/m_radius);
        parameters[5] = params(5)*(m_maxTrans/m_radius);
        transform->SetParameters(parameters);
        //std::cout<<"euler tsf parameters : "<<transform->GetParameters()<<std::endl;
        
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
            cerr<<"error while transforming mesh"<<endl;
            cerr<<e<<endl;
        }
        
        maskTsf = maskResampler->GetOutput();
        
    }
    
    return maskTsf;
}

