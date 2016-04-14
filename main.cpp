//
//  main.cpp
//  LC2_M
// Calcul de la metrique LC2
//  Created by Maxime Gérard on 17/12/15.
//  Copyright © 2015 Maxime Gérard. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>

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
#include "itkBsplineTransform.h"
#include "itkAdaptiveHistogramEqualizationImageFilter.h"

#include <cmath>

#include <dlib/optimization.h>
#include <dlib/matrix.h>

using namespace std;
using namespace dlib;




//images
typedef itk::Image<double,3> ImageType;
typedef itk::Image<unsigned char, 3> BinaryImageType;
typedef itk::Image<double,2> Image2DType;

//IO
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileReader<Image2DType> Reader2DType;
typedef itk::ImageFileReader<BinaryImageType> BinaryReaderType;
typedef itk::ImageFileWriter<ImageType> WriterType;

//iteration over images
typedef itk::ImageRegionIterator<ImageType> IteratorType;

//recalage
typedef itk::AmoebaOptimizer AmoebaOptimizerType;
typedef itk::LevenbergMarquardtOptimizer LMOptimizerType;
typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;
typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
typedef itk::Euler3DTransform<double> EulerTransformType;
typedef itk::LC2ImageToImageMetric<ImageType, ImageType> LC2MetricType;
typedef itk::TranslationTransform<double,3> TranslationType;
typedef itk::AffineTransform<double,3> AffineTransformType;
typedef itk::BSplineTransform<double,3,3> BSplineTransformType;
typedef itk::AdaptiveHistogramEqualizationImageFilter<ImageType> HistoEqualizerType;

typedef itk::LC2MVImageToImageMetric<ImageType, ImageType> LC2MVMetricType;

//pour le mapping de l'image mobile registree ac tsf determinee par registration framework
typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;

//pour la mise a l'echelle de la bande passante des intensites
typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> RescaleFilterType;
typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxCalculatorType;
typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkFilterType;

//MultiResolutionRegistration


using namespace std;


//classe pour suivre le recalage
class CommandIterationUpdate : public itk::Command
{
public :
    typedef CommandIterationUpdate Self;
    typedef itk::Command SuperClass;
    typedef itk::SmartPointer<Self> Pointer;
    itkNewMacro(Self);
    
protected:
    CommandIterationUpdate()
    {
        m_IterationNumber =0;
    }
public:
    typedef itk::AmoebaOptimizer OptimizerType;
    typedef const OptimizerType * OptimizerPointer;
    
    void Execute(itk::Object *caller, const itk::EventObject &event) ITK_OVERRIDE
    {
        Execute( (const itk::Object *)caller, event);
    }
    
    void Execute(const itk::Object * object,
                 const itk::EventObject & event) ITK_OVERRIDE
    {
        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
        if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
            return;
        }
        std::cout << m_IterationNumber++ << "   "<<endl;
        std::cout << optimizer->GetCachedValue() << "   "<<endl;
        std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
    }
private:
    unsigned long m_IterationNumber;
};




int main(int argc, const char * argv[]) {
    
    //to determine time of computation
    std::srand(time(NULL));
    std::time_t tbegin,tend;
    double texec = 0;
    tbegin = std::time(NULL);
   
    //lecture des images
    
    string filenameUS;
    string filenameIRM;
    string outputPath;
    string filenameMaskLiver;
    
    bool useLiverMask=false;
    
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
        
        if(strcmp(argv[i], "-iMaskLiver")==0)
        {
            i++;
            filenameMaskLiver= argv[i];
            useLiverMask=true;
            cout<<"Use of mask image of liver"<<endl;
        }
        
        if(strcmp(argv[i], "-o")==0)
        {
            i++;
            outputPath = argv[i];
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
    
    
    if(outputPath == "")
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
    
    
    //min max image US
    MinMaxCalculatorType::Pointer minMaxUS = MinMaxCalculatorType::New();
    minMaxUS->SetImage(image_US);
    minMaxUS->Compute();
    
    cout<<"intensity range US image : "<<"[ "<<minMaxUS->GetMinimum()<<","<<minMaxUS->GetMaximum()<<" ]"<<endl;
    
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
    
    //min max image IRM
    MinMaxCalculatorType::Pointer minMaxIRM = MinMaxCalculatorType::New();
    minMaxIRM->SetImage(image_IRM);
    minMaxIRM->Compute();
    
    cout<<"initial intensity range IRM image : "<<"[ "<<minMaxIRM->GetMinimum()<<","<<minMaxIRM->GetMaximum()<<" ]"<<endl;
    
    cout<<"done reading images"<<endl;
    
    /********************
     * Liver mask reading
     **********************/
    
    //gotta declare it outside if loop to be visible within the entire main function
    BinaryImageType::Pointer LiverMask = BinaryImageType::New();
    
    if(useLiverMask)
    {
        cout<<"Reading MRI liver mask"<<endl;
        BinaryReaderType::Pointer BinReader = BinaryReaderType::New();
        BinReader->SetImageIO(io);
        BinReader->SetFileName(filenameMaskLiver);
        try {
            BinReader->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"Error while reading liver mask"<<endl;
            cerr<<e<<endl;
        }
        
       LiverMask =BinReader->GetOutput();
    }

    
    
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

//        //verification ecriture de l'image
//            WriterType::Pointer writer6 = WriterType::New();
//            string out6 = outputPath+"/ShrunkUS.nii.gz";
//            //itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//            writer6->SetImageIO(io);
//            writer6->SetInput(US_shrunk);
//            writer6->SetFileName(out6);
//            try {
//                writer6->Update();
//            } catch (itk::ExceptionObject &e) {
//                cerr<<"error while writing rescaled image"<<endl;
//                cerr<<e<<endl;
//                return EXIT_FAILURE;
//            }
//
    cout<<"done writing shrunk US"<<std::endl;
    cout<<"taille US shrunk : "<<US_shrunk->GetLargestPossibleRegion().GetSize()<<std::endl;
    

    /*******************
     * DOWNSAMPLING IRM
     *******************/
    
    //APPROCHE MUTLIRES
    
    ShrinkFilterType::Pointer shrinkerMRI = ShrinkFilterType::New();
    shrinkerMRI->SetInput(image_IRM);
    shrinkerMRI->SetShrinkFactor(0, 2); //2 pour MRI
    shrinkerMRI->SetShrinkFactor(1, 2);
    shrinkerMRI->SetShrinkFactor(2, 2);
    try {
        shrinkerMRI->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while shrinking MRI image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    ImageType::Pointer MRI_shrunk = shrinkerMRI->GetOutput();
    
    cout<<"taille MRI shrunk : "<<MRI_shrunk->GetLargestPossibleRegion().GetSize()<<endl;
    
    /****************
     * RESCALING IRM
     ***************/
    
    //MAKE SUR MRI INTENSITIES ARE WITHIN [0,255] -> why not put the US data in the MRI range ?
    cout<<"rescaling de l'image IRM"<<endl;
    RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    //rescaler->SetInput(image_IRM);
    rescaler->SetInput(MRI_shrunk);
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
    
//    //ecrire res
//    WriterType::Pointer writer10 = WriterType::New();
//    writer10->SetImageIO(io);
//    writer10->SetInput(rescaled_IRM);
//    string out10 = outputPath+"/rescaled_MRI.nii.gz";
//    writer10->SetFileName(out10);
//    try {
//        writer10->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"Error while writing rescaled MRI"<<endl;
//        cerr<<e<<endl;
//    }

    //min max image IRM
    MinMaxCalculatorType::Pointer minMaxIRM2 = MinMaxCalculatorType::New();
    minMaxIRM2->SetImage(rescaled_IRM);
    minMaxIRM2->Compute();

    cout<<"rescaled intensity range IRM image : "<<"[ "<<minMaxIRM2->GetMinimum()<<","<<minMaxIRM2->GetMaximum()<<" ]"<<endl;
    
    
//    /*************************
//     * RECALAGE
//    *************************/
//    
//    cout<<"registration routine"<<endl;
//    
//    //test cropping MRI
//    
////    LC2_function LC2 = LC2_function(rescaled_IRM, US_shrunk);
////    LC2.limitMRI();
//    
//    //they must be column_vector !!!
//    
//    //BOBYQA initiaux FULL RES
//    //Porc 1
//    //best score : 0.183876
//    //final parameters : [-0.03514169419515513, 0.04780031905516767, 0.08510374453376125, -3.16002229212777, -0.7603111683282886, -1.5771778577960676]
//    
//    //porc 2
//    //MRI Shrunk : best score 1ere etape : 0.149085
//    //final parameters : [-0.006769331840125886, 0.09820386001984022, -0.01826162744483251, 2.7224875212753967, 0.7292929664137704, 3.0486718768199808]
//    
//    //full resolution MRI
////    best score 1ere etape : 0.156658
////    final parameters : [0.0036297748801233246, 0.013143829260491974, 0.0010727709811333918, 0.08008875517288504, -0.04607054405119398, 1.0090378953697132]
//    
//    //MRI rescaled
//    
//    //porc4
////    best score : 0.185428
////    final parameters : [-0.025424753944333096, 0.39962575090364716, 0.013544113893079517, -0.047444368532859084, -0.739287807247723, -0.5011674690804339]
//    
//    //porc6
////    best score : 0.150082
////    final parameters : [-0.0150802635531933, -0.010100620088642457, 0.004725328308930741, 8.30720384677868, 0.3776652489827957, 0.9237480980253407]
//
//
//    //MRI Shrunk + mask
////    best score 1ere etape : 0.195807
////    final parameters : [-0.004576003822187878, -0.018662112915280166, 0.021674120903660443, 1.335316780176604, 0.012001995175122276, -0.08281897027573802]
//    
//    /**********
//     *IDEAS TO TEST
//     ************/
//
//    //using mask of liver to limit ROI
//    //try with arterial phased MRA (p6)
//    //try playing on parameters to understand BOBYQA better, what about CMA-ES ?
//    //GPU programming ?
//  
//    
//    matrix<double> initialParameters (6,1);
//    initialParameters(0) = 0;
//    initialParameters(1) = 0;
//    initialParameters(2) = 0;
//    initialParameters(3) = 0;
//    initialParameters(4) = 0;
//    initialParameters(5) = 0;
//    
//    
//    //cout<<"test params intial : "<<initialParameters<<endl;
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
//    double radius = 0.9; //0.9 //zero makes no sense !!!
//    
//    //rho end
//    double precision = 0.1; //0.01 0.001
//    
//    //niter
//    double nombreIteration = 200;
//    
//    LC2_function LC2 = LC2_function(rescaled_IRM,US_shrunk,outputPath);//MRI_shrunk or rescaled
//    //make sure that the mask is computed on US_Shrunk but that we use the rescaled image to compute LC2
//    //LC2.setMovingImage(rescaled_US);
//    LC2.setMaxRot(0.3); //0.3 for best initialisation
//    LC2.setMaxTrans(10);//10
//    LC2.setRadius(radius);
//    
//    if(useLiverMask) LC2.setLiverMask(LiverMask);
//    
//   double best_score = find_max_bobyqa(LC2, initialParameters, m, x_lower, x_upper, LC2.getRadius(), precision, nombreIteration);
//    //double best_score = find_min_bobyqa(LC2_function(rescaled_IRM, US_shrunk), initialParameters, m, x_lower, x_upper, radius, precision, nombreIteration);
//    
//    cout<<"best score 1ere etape : "<<best_score<<endl;
//    //cout<<"rough parameters : "<<initialParameters<<endl;
//
////    EulerTransformType::ParametersType Step1Parameters(6);
////    
////    Step1Parameters[0] = initialParameters(0)*LC2.getMaxRot()/(LC2.getRadius());
////    Step1Parameters[1]= initialParameters(1)*LC2.getMaxRot()/(LC2.getRadius());
////    Step1Parameters[2] = initialParameters(2)*LC2.getMaxRot()/(LC2.getRadius());
////    Step1Parameters[3] = initialParameters(3)*LC2.getMaxTrans()/(LC2.getRadius());
////    Step1Parameters[4] = initialParameters(4)*LC2.getMaxTrans()/(LC2.getRadius());
////    Step1Parameters[5] = initialParameters(5)*LC2.getMaxTrans()/(LC2.getRadius());
////    cout<<"test best parameters 1ere etape: "<<Step1Parameters<<endl;
//    
//    
//    //enregistrement resultats
//    EulerTransformType::ParametersType finalParameters(6);
//    finalParameters[0] = initialParameters(0)*LC2.getMaxRot()/(LC2.getRadius());
//    finalParameters[1] = initialParameters(1)*LC2.getMaxRot()/(LC2.getRadius());
//    finalParameters[2] = initialParameters(2)*LC2.getMaxRot()/(LC2.getRadius());
//    finalParameters[3] = initialParameters(3)*LC2.getMaxTrans()/(LC2.getRadius());
//    finalParameters[4] = initialParameters(4)*LC2.getMaxTrans()/(LC2.getRadius());
//    finalParameters[5] = initialParameters(5)*LC2.getMaxTrans()/(LC2.getRadius());
//    
//    
//    //FINAL TSF
//    
//    EulerTransformType::Pointer finalTsf = EulerTransformType::New();
//    finalTsf->SetParameters(finalParameters);
//    
//    cout<<"final parameters : "<<finalTsf->GetParameters()<<endl;
//    
//    //ecritire des parametres dans un fichier txt
//    
//    string outP =outputPath+"/parameters.txt";
//    ofstream fichier(outP.c_str(),ios::out | ios::trunc );
//    
//    if(fichier)
//    {
//        fichier<<"Parameters for rigid transform : "<<endl;
//        fichier<<finalTsf->GetParameters()<<endl;
//        fichier<<" Score for this position : "<<endl;
//        fichier<<best_score<<endl;
//        fichier.close();
//    }
//    
//    else
//    {
//        cerr<<"Error in opening txt file for parameters"<<endl;
//    }
//    
//    cout<<"Writing final registered US image"<<endl;
//    
//    ImageType::SizeType sizeUS2 = US_shrunk->GetLargestPossibleRegion().GetSize();
//    ImageType::PointType origin2 = US_shrunk->GetOrigin();
//    ImageType::SpacingType spacing2 = US_shrunk->GetSpacing();
//    ImageType::PointType center2;
//    center2[0] = origin2[0]+spacing2[0]*sizeUS2[0]/2;
//    center2[1] = origin2[1]+spacing2[1]*sizeUS2[1]/2;
//    center2[2] = origin2[2]+spacing2[2]*sizeUS2[2]/2;
//    
//    
//    EulerTransformType::ParametersType eulerFixedParameters2(3);
//    eulerFixedParameters2[0] =center2[0];
//    eulerFixedParameters2[1] =center2[1];
//    eulerFixedParameters2[2] =center2[2];
//    
//    finalTsf->SetFixedParameters(eulerFixedParameters2);
//    
//    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
//    resampler->SetInput(image_US);
//    resampler->SetTransform(finalTsf);
//    resampler->SetSize(image_US->GetLargestPossibleRegion().GetSize());
//    resampler->SetOutputSpacing(image_US->GetSpacing());
//    resampler->SetOutputDirection(finalTsf->GetInverseMatrix()*image_US->GetDirection());
//    resampler->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin()));
//    //est ce que c'est une bonne idee de garder l'image full res alors qu'on fait le recalage base sur l'us downsampled ?
//    
//    cout<<"verification origine : "<<endl;
//    cout<<"avant tsf : "<<image_US->GetOrigin()<<endl;
//    cout<<"apres : "<<finalTsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin());
//    
//    ImageType::Pointer finalImage = ImageType::New();
//    finalImage = resampler->GetOutput();
//    
//    cout<<"writing final result"<<endl;
//    
//    WriterType::Pointer writer8 = WriterType::New();
//    string out8 =outputPath+"/finalregistreredUSBOBYQA.nii.gz";
//    writer8->SetImageIO(io);
//    writer8->SetInput(finalImage);
//    writer8->SetFileName(out8);
//    try {
//        writer8->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
//    
//    //target points evaluation
//    ImageType::PointType target1;
//    target1[0]=-49.000;
//    target1[1]=-23.689;
//    target1[2]=-47.019;
//    
//    ImageType::PointType target2;
//    target2[0]=-22.620;
//    target2[1]=-23.424;
//    target2[2]=-11.070;
//    
//    ImageType::PointType target3;
//    target3[0]=-44.518;
//    target3[1]=-21.495;
//    target3[2]=-50.116;
//    
//    ImageType::PointType target4;
//    target4[0]=-31.493;
//    target4[1]=1.394;
//    target4[2]=-16.948;
//    
//    ImageType::PointType target5;
//    target5[0]=-41.194;
//    target5[1]=-1.120;
//    target5[2]=-12.452;
//    
//    //transformation des targets points
//    
//    ImageType::PointType resTarget1 = finalTsf->GetInverseTransform()->TransformPoint(target1);
//    ImageType::PointType resTarget2 = finalTsf->GetInverseTransform()->TransformPoint(target2);
//    ImageType::PointType resTarget3 = finalTsf->GetInverseTransform()->TransformPoint(target3);
//    ImageType::PointType resTarget4 = finalTsf->GetInverseTransform()->TransformPoint(target4);
//    ImageType::PointType resTarget5 = finalTsf->GetInverseTransform()->TransformPoint(target5);
//    
//    //enregistrement dans fichier txt
//    
//    //ecritire des parametres dans un fichier txt
//    
//    string outT =outputPath+"/results_targets.txt";
//    ofstream fichier2(outT.c_str(),ios::out | ios::trunc );
//    
//    if(fichier)
//    {
//        fichier2<<"Transformed target points : "<<endl;
//        fichier2<<"resTarget 1 : "<<resTarget1<<endl;
//        fichier2<<"resTarget 2 : "<<resTarget2<<endl;
//        fichier2<<"resTarget 3 : "<<resTarget3<<endl;
//        fichier2<<"resTarget 4 : "<<resTarget4<<endl;
//        fichier2<<"resTarget 5 : "<<resTarget5<<endl;
//        fichier2.close();
//    }
//    
//    else
//    {
//        cerr<<"Error in opening txt file for parameters"<<endl;
//    }
    
    
    ////////////////
    //1. TRANSFORM /
    ////////////////
    
    cout<<"creation tsf"<<endl;
    
    //euler tsf
    EulerTransformType::Pointer transform = EulerTransformType::New();
    transform->SetIdentity();
    
    
    /////////////////
    // 2. OPTIMIZER /
    /////////////////
    

    //Amoeba optimizer
    
    cout<<"creation optimizer"<<endl;
    
    AmoebaOptimizerType::Pointer optimizer = AmoebaOptimizerType::New();
    
    
    //interpolateur et metrique
    
    /////////////////////
    // 3. INTERPOLATOR //
    /////////////////////
    
    cout<<"creation interpolateur"<<endl;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    
    
    /////////////////
    //  4. REGISTOR /
    ////////////////
    //pipeline de recalage
    
    RegistrationType::Pointer registration = RegistrationType::New();
    
    //settings des blocs constitutifs
   
    cout<<"settings des blocs constitutifs"<<endl;
  
    registration->SetOptimizer(optimizer);
    registration->SetTransform(transform);
    registration->SetInterpolator(interpolator);

//
//    ///////////////
//    // 4. METRIC //
//    ///////////////
    cout<<"setting metrique"<<endl;
    LC2MetricType::Pointer metric = LC2MetricType::New();

    
  //Amoeba optimizer
    
    registration->SetMetric(metric);
    
    
    //setting des images
    cout<<"setting images"<<endl;
    registration->SetFixedImage(rescaled_IRM);
    registration->SetMovingImage(US_shrunk);
    //just to be sure for the metric
    metric->SetFixed(rescaled_IRM);
    metric->ComputeGradImage();
    //metric->ComputeVesselnessImage();
    metric->SetMoving(US_shrunk);
    metric->ComputeMask();
    
    registration->SetFixedImageRegion(rescaled_IRM->GetBufferedRegion());
    
    
    RegistrationType::ParametersType initialParameters = transform->GetParameters();
    //setting des parametres optimizable
    initialParameters[0] = 0;
    initialParameters[1] = 0;
    initialParameters[2] = 0;
    initialParameters[3] = 0; //euler tsf = 6 parameters
    initialParameters[4] = 0;
    initialParameters[5] = 0;
    
    registration->SetInitialTransformParameters(initialParameters);
    
    cout<<"Initial transform param = "<<registration->GetTransform()->GetParameters()<<endl;
    
    /************
     * OPTIMIZER
     **************/
    
    //setting des params de l'optimizer
    const unsigned int numberOfParameters = transform->GetNumberOfParameters();
    AmoebaOptimizerType::ParametersType simplexDelta(numberOfParameters);
    simplexDelta[0] =0.3;
    simplexDelta[1] =0.3;
    simplexDelta[2] =0.3;
    simplexDelta[3] = 10;
    simplexDelta[4] = 10;
    simplexDelta[5] = 10;
    
    cout<<"verification simplex delta structure : "<<endl;
    cout<<simplexDelta<<endl;
    
    optimizer->AutomaticInitialSimplexOff();
    optimizer->SetInitialSimplexDelta(simplexDelta);
    
    //tolerance setting on param units and LC2 units
    optimizer->SetParametersConvergenceTolerance(0.1);
    optimizer->SetFunctionConvergenceTolerance(0.001);
    
    //max number of iterations to stop if no convergence
    optimizer->SetMaximumNumberOfIterations(200);
    
    //on le met en maximisation
    optimizer->SetMaximize(true);
    
    //Command observer
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    optimizer->AddObserver(itk::IterationEvent(), observer);
    

    cout<<"registration"<<endl;
    
        try {
            registration->Update();
            cout<<"Optimizer stop condition : "<<registration->GetOptimizer()->GetStopConditionDescription()<<endl;
        } catch (itk::ExceptionObject &e) {
            cerr<<"erreur dans le recalage"<<endl;
            cerr<<e<<endl;
            return EXIT_FAILURE;
        }
    

    RegistrationType::ParametersType finalParameters = registration->GetLastTransformParameters();
    
    double bestValue = optimizer->GetValue();
    
    //Print out results
    
    cout<<"Results : "<<endl;
    //cout<<"Translation vector : "<<"[ "<<TX<<", "<<TY<<", "<<TZ<<" ]"<<endl;
    cout<<"parametres : "<<finalParameters<<endl;
    cout<<"metric value : "<<bestValue<<endl;
    
    EulerTransformType::Pointer finalTsf = EulerTransformType::New();
    finalTsf->SetParameters(finalParameters);
    
    cout<<"final parameters : "<<finalTsf->GetParameters()<<endl;
    
    //ecritire des parametres dans un fichier txt
    
    string outP =outputPath+"/parameters.txt";
    ofstream fichier(outP.c_str(),ios::out | ios::trunc );
    
    if(fichier)
    {
        fichier<<"Parameters for rigid transform : "<<endl;
        fichier<<finalTsf->GetParameters()<<endl;
        fichier<<" Score for this position : "<<endl;
        fichier<<best_score<<endl;
        fichier.close();
    }
    
    else
    {
        cerr<<"Error in opening txt file for parameters"<<endl;
    }


    cout<<"Writing final registered US image"<<endl;
    
    ImageType::SizeType sizeUS2 = US_shrunk->GetLargestPossibleRegion().GetSize();
    ImageType::PointType origin2 = US_shrunk->GetOrigin();
    ImageType::SpacingType spacing2 = US_shrunk->GetSpacing();
    ImageType::PointType center2;
    center2[0] = origin2[0]+spacing2[0]*sizeUS2[0]/2;
    center2[1] = origin2[1]+spacing2[1]*sizeUS2[1]/2;
    center2[2] = origin2[2]+spacing2[2]*sizeUS2[2]/2;
    
    
    EulerTransformType::ParametersType eulerFixedParameters2(3);
    eulerFixedParameters2[0] =center2[0];
    eulerFixedParameters2[1] =center2[1];
    eulerFixedParameters2[2] =center2[2];
    
    finalTsf->SetFixedParameters(eulerFixedParameters2);
    
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput(image_US);
    resampler->SetTransform(finalTsf);
    resampler->SetSize(image_US->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputSpacing(image_US->GetSpacing());
    resampler->SetOutputDirection(finalTsf->GetInverseMatrix()*image_US->GetDirection());
    resampler->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin()));
    //est ce que c'est une bonne idee de garder l'image full res alors qu'on fait le recalage base sur l'us downsampled ?
    
    ImageType::Pointer finalImage = ImageType::New();
    finalImage = resampler->GetOutput();
    
    cout<<"writing final result"<<endl;
    
    WriterType::Pointer writer8 = WriterType::New();
    string out8 =outputPath+"/finalregistreredUSBOBYQA.nii.gz";
    writer8->SetImageIO(io);
    writer8->SetInput(finalImage);
    writer8->SetFileName(out8);
    try {
        writer8->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error whilte writing registered image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    
    /***********************
     * TEST TSF DEFORMABLE
     **********************/
    
    BSplineTransformType::Pointer Btransform = BSplineTransformType::New();
    
    //definition de la region de deformation
    
    
    
    
    
    std::cout << "Done running routine!\n";
  
    tend = std::time(NULL);
    texec = std::difftime(tend,tbegin);
    
    cout<<"temps en s : "<<texec<<endl;
    return 0;
}
