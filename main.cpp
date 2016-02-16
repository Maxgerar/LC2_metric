//
//  main.cpp
//  LC2_M
// Calcul de la metrique LC2
//  Created by Maxime Gérard on 17/12/15.
//  Copyright © 2015 Maxime Gérard. All rights reserved.
//

#include <iostream>
#include <string>
#include "gradient.h"
#include "LC2_function.hpp"



#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegistrationMethod.h"
#include "itkAmoebaOptimizer.h"
#include "itkNiftiImageIO.h"
#include "itkPNGImageIO.h"
#include "itkMetaImageIO.h"
#include "itkImageIterator.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRigid3Dtransform.h"
#include "itkAmoebaOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkEuler3DTransform.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include"itkLC2ImageToImageMetric.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkShrinkImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkCommand.h"
#include "itkAffineTransform.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkCastImageFilter.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkLC2MVImageToImageMetric.h"

#include <cmath>

#include <dlib/optimization.h>
#include <dlib/matrix.h>

using namespace std;
using namespace dlib;




//images
typedef itk::Image<double,3> ImageType;
typedef itk::Image<double,2> Image2DType;

//IO
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileReader<Image2DType> Reader2DType;
typedef itk::ImageFileWriter<ImageType> WriterType;

//iteration over images
typedef itk::ImageRegionIterator<ImageType> IteratorType;

//recalage
typedef itk::AmoebaOptimizer OptimizerType;
typedef itk::LevenbergMarquardtOptimizer LMOptimizerType;
typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;
typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
typedef itk::Euler3DTransform<double> EulerTransformType;
typedef itk::LC2ImageToImageMetric<ImageType, ImageType> LC2MetricType;
typedef itk::TranslationTransform<double,3> TranslationType;
typedef itk::AffineTransform<double,3> AffineTransformType;

typedef itk::LC2MVImageToImageMetric<ImageType, ImageType> LC2MVMetricType;

//pour le mapping de l'image mobile registree ac tsf determinee par registration framework
typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;

//pour la mise a l'echelle de la bande passante des intensites
typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> RescaleFilterType;
typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxCalculatorType;
typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkFilterType;

//MultiResolutionRegistration


using namespace std;

class CommandIterationUpdate : public itk::Command
{
public:
    typedef CommandIterationUpdate  Self;
    typedef itk::Command            Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    itkNewMacro( Self );
protected:
    CommandIterationUpdate() {};
public:
    typedef itk::LevenbergMarquardtOptimizer     OptimizerType;
    typedef const OptimizerType *                OptimizerPointer;
    void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
    {
        Execute( (const itk::Object *)caller, event);
    }
    void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
    {
        OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( object );
        if( optimizer == ITK_NULLPTR )
        {
            itkExceptionMacro( "Could not cast optimizer." );
        }
        if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
            return;
        }
        //std::cout << "Value = " << optimizer->GetCachedValue() << std::endl;
        std::cout << "Position = "  << optimizer->GetCachedCurrentPosition();
        std::cout << std::endl << std::endl;
    }
};

////classe pour suivre le recalage
//class CommandIterationUpdate : public itk::Command
//{
//public :
//    typedef CommandIterationUpdate Self;
//    typedef itk::Command SuperClass;
//    typedef itk::SmartPointer<Self> Pointer;
//    itkNewMacro(Self);
//    
//protected:
//    CommandIterationUpdate()
//    {
//        m_IterationNumber =0;
//    }
//public:
//    typedef itk::AmoebaOptimizer OptimizerType;
//    typedef const OptimizerType * OptimizerPointer;
//    
//    void Execute(itk::Object *caller, const itk::EventObject &event) ITK_OVERRIDE
//    {
//        Execute( (const itk::Object *)caller, event);
//    }
//    
//    void Execute(const itk::Object * object,
//                 const itk::EventObject & event) ITK_OVERRIDE
//    {
//        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
//        if( ! itk::IterationEvent().CheckEvent( &event ) )
//        {
//            return;
//        }
//        std::cout << m_IterationNumber++ << "   "<<endl;
//        std::cout << optimizer->GetCachedValue() << "   "<<endl;
//        std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
//    }
//private:
//    unsigned long m_IterationNumber;
//};

////classe pour suivre le recalage multires
//template <typename TRegistration>
//class RegistrationInterfaceCommand : public itk::Command
//{
//public:
//    typedef  RegistrationInterfaceCommand   Self;
//    typedef  itk::Command                   Superclass;
//    typedef  itk::SmartPointer<Self>        Pointer;
//    itkNewMacro( Self );
//protected:
//    RegistrationInterfaceCommand() {};
//public:
//    typedef   TRegistration                              RegistrationType;
//    typedef   RegistrationType *                         RegistrationPointer;
//    typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType;
//    typedef   OptimizerType *                            OptimizerPointer;
//    void Execute(itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
//    {
//        if( !(itk::IterationEvent().CheckEvent( &event )) )
//        {
//            return;
//        }
//        RegistrationPointer registration = static_cast<RegistrationPointer>( object );
//        if(registration == ITK_NULLPTR)
//        {
//            return;
//        }
//        OptimizerPointer optimizer = static_cast< OptimizerPointer >(registration->GetModifiableOptimizer() );
//        std::cout << "-------------------------------------" << std::endl;
//        std::cout << "MultiResolution Level : "
//        << registration->GetCurrentLevel()  << std::endl;
//        std::cout << std::endl;
//        if ( registration->GetCurrentLevel() == 0 )
//        {
//            optimizer->SetMaximumStepLength( 16.00 );
//            optimizer->SetMinimumStepLength( 0.01 );
//        }
//        else
//        {
//            optimizer->SetMaximumStepLength( optimizer->GetMaximumStepLength() / 4.0 );
//            optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength() / 10.0 );
//        }
//    }
//    void Execute(const itk::Object * , const itk::EventObject & ) ITK_OVERRIDE
//    { return; }
//};
//
//class CommandIterationUpdate : public itk::Command
//{
//public:
//    typedef  CommandIterationUpdate   Self;
//    typedef  itk::Command             Superclass;
//    typedef  itk::SmartPointer<Self>  Pointer;
//    itkNewMacro( Self );
//protected:
//    CommandIterationUpdate() {};
//public:
//    typedef   itk::RegularStepGradientDescentOptimizer OptimizerType;
//    typedef   const OptimizerType *                    OptimizerPointer;
//    void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
//    {
//        Execute( (const itk::Object *)caller, event);
//    }
//    void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
//    {
//        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
//        if( !(itk::IterationEvent().CheckEvent( &event )) )
//        {
//            return;
//        }
//        std::cout << optimizer->GetCurrentIteration() << "   ";
//        std::cout << optimizer->GetValue() << "   ";
//        std::cout << optimizer->GetCurrentPosition() << std::endl;
//    }
//};

int main(int argc, const char * argv[]) {
    
    //to determine time of computation
    std::srand(time(NULL));
    std::time_t tbegin,tend;
    double texec = 0;
    tbegin = std::time(NULL);
   
    //lecture des images
    
    string filenameUS;
    string filenameIRM;
    string outputFilename;
    
    /*********************
     * RUNTIME ARG
     ******************/
    
    cout<<"runtime arguments acquisition"<<endl;
    
    for(int i = 0; i < argc; ++i)
    {
        //input image US
        if(strcmp(argv[i], "-iUS")==0)
        {
            i++;
            filenameUS = argv[i];
        }
        
        //input image IRM
        if(strcmp(argv[i], "-iIRM")==0)
        {
            i++;
            filenameIRM = argv[i];
        }
        
        if(strcmp(argv[i], "-o")==0)
        {
            i++;
            outputFilename = argv[i];
        }
        
        
    }
    
    /**********************
     * VERIFICATION INPUTS
     ********************/
    
    cout<<"input error handling"<<endl;
    
    if(filenameUS == "")
    {
        cerr<<"Input US file not provided"<<endl;
        return EXIT_FAILURE;
        
    }
    
    if(filenameIRM == "")
    {
        cerr<<"Input MRI file not provided"<<endl;
        return EXIT_FAILURE;
    }
    
    if(outputFilename == "")
    {
        cerr<<"output path not provided"<<endl;
        return EXIT_FAILURE;
    }
    
    /**********************
     * US READING
     *********************/
    
    cout<<"Reading images"<<endl;
    
    ImageType::Pointer image_US = ImageType::New();
    ImageType::Pointer image_IRM = ImageType::New();
    
    
    ReaderType::Pointer reader1 = ReaderType::New();
    itk::MetaImageIO::Pointer m_io = itk::MetaImageIO::New();
    reader1->SetImageIO(m_io);
    reader1->SetFileName(filenameUS);
    try {
        reader1->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"Error while reading US image"<<endl;
        cerr<<e<<endl;
        EXIT_FAILURE;
    }
    
    image_US = reader1->GetOutput();
    cout<<"test lecture US"<<endl;
    cout<<"dimensions US : "<<image_US->GetLargestPossibleRegion().GetSize()<<endl;
    
    
//    //min max image US
//    MinMaxCalculatorType::Pointer minMaxUS = MinMaxCalculatorType::New();
//    minMaxUS->SetImage(image_US);
//    minMaxUS->Compute();
//    
//    cout<<"intensity range US image : "<<"[ "<<minMaxUS->GetMinimum()<<","<<minMaxUS->GetMaximum()<<" ]"<<endl;
    
    /**********************
     * IRM READING
     *********************/
    
    ReaderType::Pointer reader2 = ReaderType::New();
    itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
    reader2->SetImageIO(io);
    reader2->SetFileName(filenameIRM);
    try {
        reader2->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"Error while reading IRM image"<<endl;
        cerr<<e<<endl;
        EXIT_FAILURE;
    }
    
    image_IRM = reader2->GetOutput();
    cout<<"test lecture IRM"<<endl;
    cout<<"dimensions IRM : "<<image_IRM->GetLargestPossibleRegion().GetSize()<<endl;
    
//    //min max image IRM
//    MinMaxCalculatorType::Pointer minMaxIRM = MinMaxCalculatorType::New();
//    minMaxIRM->SetImage(image_IRM);
//    minMaxIRM->Compute();
//    
//    cout<<"initial intensity range IRM image : "<<"[ "<<minMaxIRM->GetMinimum()<<","<<minMaxIRM->GetMaximum()<<" ]"<<endl;
    
    cout<<"done reading images"<<endl;
    
    /*******************
     * DOWNSAMPLING US
     ******************/
    //RESOLUTION ADAPTATION US
    ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
    shrinkFilter->SetInput(image_US);
    
    //recuperation des spacing des deux images pour savoir quel facteur utiliser
    ImageType::SpacingType spacingUS = image_US->GetSpacing();
    ImageType::SpacingType spacingIRM = image_IRM->GetSpacing();
    //calcul du spacing factor
    int shrinkX = int(2*spacingIRM[0]/spacingUS[0]);
    int shrinkY = int(2*spacingIRM[1]/spacingUS[1]);
    int shrinkZ = int(2*spacingIRM[2]/spacingUS[2]);
    
    cout<<"shrinking factors : "<<shrinkX<<", "<<shrinkY<<", "<<shrinkZ<<endl;

    //porc 6 : (3,3,2)/ porc 1 (5,4,3)
    shrinkFilter->SetShrinkFactor(0, shrinkX);
    shrinkFilter->SetShrinkFactor(1, shrinkY);
    shrinkFilter->SetShrinkFactor(2, shrinkZ);
    try {
        shrinkFilter->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while downsampling US image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    ImageType::Pointer US_shrunk = shrinkFilter->GetOutput();
//
        //verification ecriture de l'image
            WriterType::Pointer writer6 = WriterType::New();
            string out6 = "/Users/maximegerard/Documents/ShrunkUS.nii.gz";
          //itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
            writer6->SetImageIO(io);
            writer6->SetInput(US_shrunk);
            writer6->SetFileName(out6);
            try {
                writer6->Update();
            } catch (itk::ExceptionObject &e) {
                cerr<<"error while writing rescaled image"<<endl;
                cerr<<e<<endl;
                return EXIT_FAILURE;
            }
    
    cout<<"done writing shrunk US"<<endl;
    
    /****************
    * RESCALING IRM
     ***************/
    
    cout<<"rescaling de l'image IRM"<<endl;
    RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetInput(image_IRM);
    rescaler->SetOutputMinimum(0);
    rescaler->SetOutputMaximum(255);
    try {
        rescaler->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while rescaling image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    
    ImageType::Pointer rescaled_IRM = rescaler->GetOutput();
    
//    //min max image IRM
//    MinMaxCalculatorType::Pointer minMaxIRM2 = MinMaxCalculatorType::New();
//    minMaxIRM2->SetImage(rescaled_IRM);
//    minMaxIRM2->Compute();
//    
//    cout<<"rescaled intensity range IRM image : "<<"[ "<<minMaxIRM2->GetMinimum()<<","<<minMaxIRM2->GetMaximum()<<" ]"<<endl;
    
//    WriterType::Pointer writer6 = WriterType::New();
//    string out6 = "/Users/maximegerard/Documents/IRMRescaled.nii.gz";
//    writer6->SetImageIO(io);
//    writer6->SetInput(rescaled_IRM);
//    writer6->SetFileName(out6);
//    try {
//        writer6->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error while writing rescaled image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
//   
    
    
    /*************************
     * RECALAGE
    *************************/
    
    cout<<"registration routine"<<endl;
    
    //test autre registration
    
    //they must be column_vector !!!
    
//    matrix<double> initialParameters (6,1);
//    initialParameters(0) = 0;
//    initialParameters(1) = 0;
//    initialParameters(2) = 0;
//    initialParameters(3) = 0;
//    initialParameters(4) = 0;
//    initialParameters(5) = 0;
//    
//    cout<<"test params intial : "<<initialParameters<<endl;
//    
//    matrix<double> x_lower (6,1); //-0.5 rot, -10 trans
//    x_lower(0) = -1;
//    x_lower(1) = -1;
//    x_lower(2) = -1;
//    x_lower(3) = -1;
//    x_lower(4) = -1;
//    x_lower(5) = -1;
//    
//    matrix<double> x_upper (6,1); //0.5 rot, 10 trans
//    x_upper(0) = 1;
//    x_upper(1) = 1;
//    x_upper(2) = 1;
//    x_upper(3) = 1;
//    x_upper(4) = 1;
//    x_upper(5) = 1;
//    
//    double m = 13;
//    
//    //rho begin
//    double radius = 0.1; //zero makes no sense !!!
//    
//    //rho end
//    double precision = 0.001;
//    
//    //niter
//    double nombreIteration = 200;
//    
//   double best_score = find_max_bobyqa(LC2_function(rescaled_IRM,US_shrunk), initialParameters, m, x_lower, x_upper, radius, precision, nombreIteration);
//    //double best_score = find_min_bobyqa(LC2_function(rescaled_IRM, US_shrunk), initialParameters, m, x_lower, x_upper, radius, precision, nombreIteration);
//    
//    cout<<"best score : "<<best_score<<endl;
//    cout<<"test best parameters : "<<initialParameters<<endl; // ok faut juste pas oublier de tout mutliplier par 5 pour les rotation et 100 pour les translations
    
    
    
//    
//    ////////////////
//    //1. TRANSFORM /
//    ////////////////
//    
//    cout<<"creation tsf"<<endl;
//    
//    //euler tsf
//    EulerTransformType::Pointer transform = EulerTransformType::New();
//    transform->SetIdentity();
//    
//    //tsl transform
//    //TranslationType::Pointer transform = TranslationType::New();
//    
//    
//    /////////////////
//    // 2. OPTIMIZER /
//    /////////////////
//    
// //LM optimizer
//    cout<<"creation LM optimizer"<<endl;
//    LMOptimizerType::Pointer optimizer = LMOptimizerType::New();
//    optimizer->SetUseCostFunctionGradient(false);
//    
//    //Amoeba optimizer
//    
////    cout<<"creation optimizer"<<endl;
////    
////    OptimizerType::Pointer optimizer = OptimizerType::New();
//    
//    
//    //interpolateur et metrique
//    
//    /////////////////////
//    // 3. INTERPOLATOR //
//    /////////////////////
//    
////    cout<<"creation interpolateur"<<endl;
////    InterpolatorType::Pointer interpolator = InterpolatorType::New();
//    
//    
//    /////////////////
//    //  4. REGISTOR /
//    ////////////////
//    //pipeline de recalage
//    
//    //RegistrationType::Pointer registration = RegistrationType::New();
//    
//    //settings des blocs constitutifs
//   
//    cout<<"settings des blocs constitutifs"<<endl;
//  
////    registration->SetOptimizer(optimizer);
////    registration->SetTransform(transform);
////    registration->SetInterpolator(interpolator);
//
//    
//    ///////////////
//    // 4. METRIC //
//    ///////////////
//    cout<<"setting metrique"<<endl;
//    //LC2MetricType::Pointer metric = LC2MetricType::New();
//    
//    //test metric mutliple valued
//    LC2MVMetricType::Pointer metricMV = LC2MVMetricType::New();
//    metricMV->SetMoving(US_shrunk);
//    metricMV->SetFixed(rescaled_IRM);
//    metricMV->ComputeGradImage();
//    metricMV->ComputeMask();
//    
//    EulerTransformType::ParametersType parameters(6);
//    parameters[0]=0;
//    parameters[1]=0;
//    parameters[2]=0;
//    parameters[3]=0;
//    parameters[4]=0;
//    parameters[5]=0;
//    transform->SetParameters(parameters);
//    metricMV->SetTransform(transform);
//    cout<<"params initiaux : "<<metricMV->GetTransform()->GetParameters()<<endl;
//    metricMV->SetTransformParameters(parameters);
//    //metricMV->GetValue(parameters);
//    
//    
//    OptimizerType::ScalesType scales(transform->GetNumberOfParameters());
//    //setting the scales and ranges for the parameters of the transform
//    const double translationScale = 20.0; //dynamic range for tsl
//    const double rotationScale = 1.0; //dynamic range for rotations
//    
//    scales[0] = 1/rotationScale;
//    scales[1] = 1/rotationScale;
//    scales[2] = 1/rotationScale;
//    scales[3] = 1/translationScale;
//    scales[4] = 1/translationScale;
//    scales[5] = 1/translationScale;
//    
//    unsigned long numberOfIterations = 200;
//    //double gradientTolerance = 1e-4;
//    double valueTolerance = 1e-3;
//    double epsilonFunction = 1e-4;
//    
//    optimizer->SetScales(scales);
//    optimizer->SetNumberOfIterations(numberOfIterations);
//    optimizer->SetValueTolerance(valueTolerance);
//    //optimizer->SetGradientTolerance(gradientTolerance);
//    optimizer->SetEpsilonFunction(epsilonFunction);
//    optimizer->SetCostFunction(metricMV);
//    optimizer->SetInitialPosition(parameters);
//    
//    //optimizer->UseCostFunctionGradientOff();
//    
//    //observateur
//    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
//    optimizer->AddObserver(itk::IterationEvent(), observer);
//    
//    cout<<"start optimisation "<<endl;
//    //optimizer->StartOptimization();
//    
//    try {
//        optimizer->StartOptimization();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error while optimizing cost function"<<endl;
//        cerr<<e<<endl;
//    }
//    
//    
//    cout<<"transform best param according to optimizer : "<<optimizer->GetCurrentPosition()<<endl;
//    cout<<"best param test transform : "<<transform->GetParameters()<<endl;
    //cout<<"best value : "<<optimizer->GetValue()<<endl;


    
    
//    registration->SetMetric(metric);
//    
//    
//    //setting des images
//    cout<<"setting images"<<endl;
//    registration->SetFixedImage(rescaled_IRM);
//    registration->SetMovingImage(US_shrunk);
//    //just to be sure for the metric
//    metric->SetFixed(rescaled_IRM);
//    metric->ComputeGradImage();
//    //metric->ComputeVesselnessImage();
//    metric->SetMoving(US_shrunk);
//    metric->ComputeMask();
//    
//    registration->SetFixedImageRegion(rescaled_IRM->GetBufferedRegion());
//    
//    
//    RegistrationType::ParametersType initialParameters = transform->GetParameters();
//    //setting des parametres optimizable
//    initialParameters[0] = 0.0;
//    initialParameters[1] = 0.0;
//    initialParameters[2] = 0.0;
//    initialParameters[3] = 0.0; //euler tsf = 6 parameters
//    initialParameters[4] = 0.0;
//    initialParameters[5] = 0.0;
//    
//    registration->SetInitialTransformParameters(initialParameters);
//    
//    cout<<"Initial transform param = "<<registration->GetTransform()->GetParameters()<<endl;
//    
//    /************
//     * OPTIMIZER
//     **************/
//    
//    //setting des params de l'optimizer
//    const unsigned int numberOfParameters = transform->GetNumberOfParameters();
//    OptimizerType::ParametersType simplexDelta(numberOfParameters);
//    simplexDelta[0] =0.5;
//    simplexDelta[1] =0.5;
//    simplexDelta[2] =0.5;
//    simplexDelta[3] = 10;
//    simplexDelta[4] = 10;
//    simplexDelta[5] = 10;
//    
//    cout<<"verification simplex delta structure : "<<endl;
//    cout<<simplexDelta<<endl;
//    
//    optimizer->AutomaticInitialSimplexOff();
//    optimizer->SetInitialSimplexDelta(simplexDelta);
//    
//    //tolerance setting on param units and LC2 units
//    optimizer->SetParametersConvergenceTolerance(0.1);
//    optimizer->SetFunctionConvergenceTolerance(0.001);
//    
//    //max number of iterations to stop if no convergence
//    optimizer->SetMaximumNumberOfIterations(200);
//    
//    //on le met en maximisation
//    optimizer->SetMaximize(true);
//    
//    //Command observer
//    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
//    optimizer->AddObserver(itk::IterationEvent(), observer);
//    
//
//    cout<<"registration"<<endl;
//    
//        try {
//            registration->Update();
//            cout<<"Optimizer stop condition : "<<registration->GetOptimizer()->GetStopConditionDescription()<<endl;
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"erreur dans le recalage"<<endl;
//            cerr<<e<<endl;
//            return EXIT_FAILURE;
//        }
//    
//
//    RegistrationType::ParametersType finalParameters = registration->GetLastTransformParameters();
//    
//    double bestValue = optimizer->GetValue();
//    
//    //Print out results
//    
//    cout<<"Results : "<<endl;
//    //cout<<"Translation vector : "<<"[ "<<TX<<", "<<TY<<", "<<TZ<<" ]"<<endl;
//    cout<<"parametres : "<<finalParameters<<endl;
//    cout<<"metric value : "<<bestValue<<endl;
//    
//    cout<<"formation de l'image recalee"<<endl;
//    
//    //TranslationType::Pointer finalTsl = TranslationType::New();
//    TransformType::Pointer finalTsf = TransformType::New();
//    
//    finalTsf->SetParameters(finalParameters);
//    finalTsf->SetFixedParameters(registration->GetTransform()->GetFixedParameters());
//    
//    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
//    resampler->SetTransform(finalTsf);
//    resampler->SetInput(image_US);
//    resampler->SetOutputDirection(finalTsf->GetInverseMatrix()*image_US->GetDirection());
//    resampler->SetSize(image_US->GetLargestPossibleRegion().GetSize());
//    resampler->SetOutputSpacing(image_US->GetSpacing());
//    resampler->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin()));
//  
//    ImageType::Pointer outputImage = ImageType::New();
//    outputImage = resampler->GetOutput();
//    
//    WriterType::Pointer writer7 = WriterType::New();
//    string out7 ="/Users/maximegerard/Documents/finalregistreredUS.nii.gz";
//    writer7->SetImageIO(io);
//    writer7->SetInput(outputImage);
//    writer7->SetFileName(out7);
//    try {
//        writer7->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
    
    
    
    
    
    std::cout << "Done running routine!\n";
  
    tend = std::time(NULL);
    texec = std::difftime(tend,tbegin);
    
    cout<<"temps en s : "<<texec<<endl;
    return 0;
}
