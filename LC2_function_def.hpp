//
//  LC2_function_def.hpp
//  LC2_M
//
//  Created by Maxime Gérard on 26/03/16.
//  Copyright © 2016 Maxime Gérard. All rights reserved.
//

#ifndef LC2_function_def_hpp
#define LC2_function_def_hpp

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
//#include "itkEuler3DTransform.h"
#include "itkBsplineTransform.h"
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
//typedef itk::Euler3DTransform<double> EulerTransformType;
typedef itk::BSplineTransform<double,3,3> BsplineTsfType;
typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
typedef itk::ResampleImageFilter<MaskType, MaskType> ResamplerBinaryType;
typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkFilterType;
typedef itk::ShrinkImageFilter<MaskType, MaskType> BinaryShrinkFilterType;

#endif /* LC2_function_def_hpp */
