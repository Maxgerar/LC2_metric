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
#ifndef itkLC2ImageToImageMetric_h
#define itkLC2ImageToImageMetric_h

#include "itkImageToImageMetric.h"
#include "itkimage.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "itkNeighborhoodIterator.h"
#include "gradient.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkNiftiImageIO.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkImageRegionIterator.h"
#include "itkExtractImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include <time.h>
#include <vector>
#include "itkResampleImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkEuler3DTransform.h"
#include "itkShrinkImageFilter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"




namespace itk
{
template< typename TFixedImage, typename TMovingImage >
class LC2ImageToImageMetric:
  public ImageToImageMetric< TFixedImage, TMovingImage >
{
public:

  /** Standard class typedefs. */
  typedef LC2ImageToImageMetric                           Self;
  typedef ImageToImageMetric< TFixedImage, TMovingImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LC2ImageToImageMetric, Object);

  /** Types transferred from the base class */
  typedef typename Superclass::RealType                RealType;
  typedef typename Superclass::TransformType           TransformType;
  typedef typename Superclass::TransformPointer        TransformPointer;
  typedef typename Superclass::TransformParametersType TransformParametersType;
  typedef typename Superclass::TransformJacobianType   TransformJacobianType;
  typedef typename Superclass::GradientPixelType       GradientPixelType;
  typedef typename Superclass::OutputPointType         OutputPointType;
  typedef typename Superclass::InputPointType          InputPointType;

  typedef typename Superclass::MeasureType             MeasureType;
  typedef typename Superclass::DerivativeType          DerivativeType;
  typedef typename Superclass::FixedImageType          FixedImageType;
  typedef typename Superclass::MovingImageType         MovingImageType;
  typedef typename Superclass::FixedImageConstPointer  FixedImageConstPointer;
  typedef typename Superclass::MovingImageConstPointer MovingImageConstPointer;
    

    //matrices pour la metrique
  typedef Matrix<double, 343,3> MType;
  typedef Matrix<double,3,343> MTType;
  typedef Matrix<double,3,3> M3Type;
  typedef Vector<double, 343> VType;
  typedef Vector<double,3> PType;
    
    //for gradient computation
    typedef GradientMagnitudeImageFilter<TFixedImage, TFixedImage> GradientFilterType;
    typedef ImageFileWriter<TFixedImage> WriterType;
    
    //mask for MRI
    typedef Image<unsigned char, 3> MaskType;
    typedef ImageFileWriter<MaskType> BinaryWriterType;
    typedef BinaryThresholdImageFilter<TMovingImage, MaskType> BinaryThresholdFilterType;
 
    
    typedef FlatStructuringElement<3> kernelType;
    typedef BinaryMorphologicalClosingImageFilter<MaskType, MaskType, kernelType> CloserType;
    typedef BinaryBallStructuringElement<unsigned char,3> StructuringElementType;
    typedef BinaryDilateImageFilter<MaskType, MaskType, StructuringElementType> DilataterType;

    //cropping
    typedef ExtractImageFilter<TFixedImage, TFixedImage> ExtractorType;
    typedef ExtractImageFilter<MaskType, MaskType> BinaryExtractorType;
 
    
    //image iterator
  
    typedef NeighborhoodIterator<TMovingImage> NeighborhoodIteratorType;
    typedef ImageRegionConstIterator<TMovingImage> ImageConstIteratorType;
    typedef ImageRegionIterator<TMovingImage> ImageIteratorType;
    typedef ImageRegionIterator<MaskType> BinaryImageIteratorType;
    
    //interpolator
    typedef LinearInterpolateImageFunction<TMovingImage,double> LinearInterpolatorType;
    
    //vesselness filtering
    typedef itk::SymmetricSecondRankTensor<double,3> HessianPixelType;
    typedef itk::Image<HessianPixelType,3> HessianImageType;
    typedef itk::Hessian3DToVesselnessMeasureImageFilter<double> VesselnessMeasureFilterType;
    typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType, HessianImageType,ImageType> MultiScaleEnhancementFilterType;
    
    //for image tsf
    typedef TranslationTransform<double,3> TranslationTransformType;
    typedef Euler3DTransform<double> EulerTransformType;
    typedef ResampleImageFilter<TMovingImage, TMovingImage> ResamplerType;
    typedef ResampleImageFilter<MaskType, MaskType> ResamplerBinaryType;
    typedef ShrinkImageFilter<TMovingImage, TMovingImage> ShrinkFilterType;
    typedef ShrinkImageFilter<MaskType, MaskType> BinaryShrinkFilterType;

  /** Get the derivatives of the match measure. */
  void GetDerivative(const TransformParametersType & parameters,
                     DerivativeType & Derivative) const ITK_OVERRIDE;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue(const TransformParametersType & parameters) const ITK_OVERRIDE;
    //const TransformParametersType & parameters

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative(const TransformParametersType & parameters,
                             MeasureType & Value, DerivativeType & Derivative) const ITK_OVERRIDE;
  
    double test(const TransformParametersType &parameters);
    

  void ComputeMask();

  /** Set/Get SubtractMean boolean. If true, the sample mean is subtracted
   * from the sample values in the cross-correlation formula and
   * typically results in narrower valleys in the cost function.
   * Default value is false. */
  itkSetMacro(lc2final, double);
  itkGetConstReferenceMacro(lc2final, double);
  //itkBooleanMacro(SubtractMean);
    
  void SetImages(typename TFixedImage::Pointer US,typename TMovingImage::Pointer IRM);
  void SetFixed(typename TFixedImage::Pointer IRM);
  void SetMoving(typename TMovingImage::Pointer US);
    
    void ComputeGradImage();
    void ComputeVesselnessImage();


protected:
  LC2ImageToImageMetric();
  virtual ~LC2ImageToImageMetric() {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
    LC2ImageToImageMetric(const Self &); //purposely not
                                                         // implemented
    void operator=(const Self &);                         //purposely not
                                                         // implemented

  //bool m_SubtractMean;
  double m_lc2final;
  double variancesum;
  double lc2varsum;
  MaskType::Pointer m_mask;
  typename TFixedImage::Pointer m_grad;




};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLC2ImageToImageMetric.hxx"
#endif

#endif
