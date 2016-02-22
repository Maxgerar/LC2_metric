//
//  LC2_function.cpp
//  LC2_M
//
//  Created by Maxime Gérard on 10/02/16.
//  Copyright © 2016 Maxime Gérard. All rights reserved.
//

#include "LC2_function.hpp"

LC2_function::LC2_function(ImageType::Pointer im_Fixed, ImageType::Pointer im_Moving)
{
    std::cout<<"creating measure object"<<std::endl;
    m_FixedImage = im_Fixed;
    m_MovingImage = im_Moving;
    computeGradient();
    computeMask();
    track =0;
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
//    std::string out3 = "/Users/maximegerard/Documents/testmask.nii.gz";
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
    std::cout<<"done writing final mask image"<<std::endl;

    
}


