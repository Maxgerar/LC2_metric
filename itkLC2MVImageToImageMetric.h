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
#ifndef itkLC2MVImageToImageMetric_h
#define itkLC2MVImageToImageMetric_h

//#include "itkImageToImageMetric.h"
#include "itkMultipleValuedCostFunction.h"
#include "itkTransform.h"
#include "itkPoint.h"
#include "itkImageToImageMetric.h"
#include "itkimage.h"
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
#include "itkArray.h"



namespace itk
{
/** \class LC2ImageToImageMetric
 * \brief Computes similarity between two images to be registered
 *
  */
template< typename TFixedImage, typename TMovingImage >
class LC2MVImageToImageMetric:
  public MultipleValuedCostFunction
{
public:

  /** Standard class typedefs. */
  typedef LC2MVImageToImageMetric         Self;
  typedef MultipleValuedCostFunction   Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;
  
    typedef Superclass::ParametersValueType CoordinateRepresentationType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LC2MVImageToImageMetric, Object);

  /** Types transferred from the base class */
//  typedef typename Superclass::RealType                RealType;
//  typedef typename Superclass::TransformType           TransformType;
//  typedef typename Superclass::TransformPointer        TransformPointer;
//  typedef typename Superclass::TransformParametersType TransformParametersType;
//  typedef typename Superclass::TransformJacobianType   TransformJacobianType;
//  typedef typename Superclass::GradientPixelType       GradientPixelType;
//  typedef typename Superclass::OutputPointType         OutputPointType;
//  typedef typename Superclass::InputPointType          InputPointType;
    
    typedef Transform<CoordinateRepresentationType,3,3> TransformType;
    typedef typename TransformType::Pointer TransformPointer;
    typedef typename TransformType::InputPointType InputPointType;
    typedef typename TransformType::OutputPointType OutputPointType;
    typedef typename TransformType::ParametersType TransformParametersType;
    typedef typename TransformType::JacobianType TransformJacobianType;

  typedef typename Superclass::MeasureType             MeasureType;
  typedef typename Superclass::DerivativeType          DerivativeType;
    typedef Superclass::ParametersType  ParametersType;
    
    itkSetObjectMacro(Transform,TransformType);
    itkGetModifiableObjectMacro(Transform,TransformType);
    
//  typedef typename Superclass::FixedImageType          FixedImageType;
//  typedef typename Superclass::MovingImageType         MovingImageType;
//  typedef typename Superclass::FixedImageConstPointer  FixedImageConstPointer;
//  typedef typename Superclass::MovingImageConstPointer MovingImageConstPointer;
    
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
    
    //image tsf
    typedef Euler3DTransform<double> EulerTransformType;
    typedef ResampleImageFilter<TMovingImage, TMovingImage> ResamplerType;
    typedef ResampleImageFilter<MaskType, MaskType> ResamplerBinaryType;
    typedef ShrinkImageFilter<TMovingImage, TMovingImage> ShrinkFilterType;
    typedef ShrinkImageFilter<MaskType, MaskType> BinaryShrinkFilterType;

//    //return value
//    itk::Array<double> ArrayType;

  /** Get the derivatives of the match measure. */
  void GetDerivative(const TransformParametersType & parameters,
                     DerivativeType & Derivative) const ITK_OVERRIDE;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue(const TransformParametersType & parameters) const ITK_OVERRIDE;
    

  unsigned int GetNumberOfValues() const ITK_OVERRIDE;
  void SetNumberOfValues(const int number) const;

    
    void ComputeMask();
    void ComputeGradImage();
    void SetImages(typename TFixedImage::Pointer US,typename TMovingImage::Pointer IRM);
    void SetFixed(typename TFixedImage::Pointer IRM);
    void SetMoving(typename TMovingImage::Pointer US);
    
    void SetTransformParameters(const ParametersType &parameters) const;
    void SetTransform(const TransformPointer &transform) const;
    
    virtual unsigned int GetNumberOfParameters(void) const ITK_OVERRIDE
    {
        return m_Transform->GetNumberOfParameters();
    }


  /** Set/Get SubtractMean boolean. If true, the sample mean is subtracted
   * from the sample values in the cross-correlation formula and
   * typically results in narrower valleys in the cost function.
   * Default value is false. */
//  itkSetMacro(SubtractMean, bool);
//  itkGetConstReferenceMacro(SubtractMean, bool);
//  itkBooleanMacro(SubtractMean);

protected:
  LC2MVImageToImageMetric();
  virtual ~LC2MVImageToImageMetric() {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;
    
  mutable TransformPointer m_Transform;


private:
  LC2MVImageToImageMetric(const Self &); //purposely not
                                                         // implemented
  void operator=(const Self &);                          //purposely not
                                                         // implemented

    MaskType::Pointer m_mask;
    unsigned int m_numberofValues;
    typename TFixedImage::Pointer m_grad;
    typename TFixedImage::Pointer m_fixedImage;
    typename TMovingImage::Pointer m_movingImage;
  
  //bool m_SubtractMean;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLC2MVImageToImageMetric.hxx"
#endif

#endif
