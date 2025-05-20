// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME particleDic
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "particle.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_particle(void *p = nullptr);
   static void *newArray_particle(Long_t size, void *p);
   static void delete_particle(void *p);
   static void deleteArray_particle(void *p);
   static void destruct_particle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::particle*)
   {
      ::particle *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::particle >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("particle", ::particle::Class_Version(), "particle.h", 8,
                  typeid(::particle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::particle::Dictionary, isa_proxy, 4,
                  sizeof(::particle) );
      instance.SetNew(&new_particle);
      instance.SetNewArray(&newArray_particle);
      instance.SetDelete(&delete_particle);
      instance.SetDeleteArray(&deleteArray_particle);
      instance.SetDestructor(&destruct_particle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::particle*)
   {
      return GenerateInitInstanceLocal((::particle*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::particle*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr particle::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *particle::Class_Name()
{
   return "particle";
}

//______________________________________________________________________________
const char *particle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::particle*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int particle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::particle*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *particle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::particle*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *particle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::particle*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void particle::Streamer(TBuffer &R__b)
{
   // Stream an object of class particle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(particle::Class(),this);
   } else {
      R__b.WriteClassBuffer(particle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_particle(void *p) {
      return  p ? new(p) ::particle : new ::particle;
   }
   static void *newArray_particle(Long_t nElements, void *p) {
      return p ? new(p) ::particle[nElements] : new ::particle[nElements];
   }
   // Wrapper around operator delete
   static void delete_particle(void *p) {
      delete ((::particle*)p);
   }
   static void deleteArray_particle(void *p) {
      delete [] ((::particle*)p);
   }
   static void destruct_particle(void *p) {
      typedef ::particle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::particle

namespace ROOT {
   static TClass *vectorlEparticlegR_Dictionary();
   static void vectorlEparticlegR_TClassManip(TClass*);
   static void *new_vectorlEparticlegR(void *p = nullptr);
   static void *newArray_vectorlEparticlegR(Long_t size, void *p);
   static void delete_vectorlEparticlegR(void *p);
   static void deleteArray_vectorlEparticlegR(void *p);
   static void destruct_vectorlEparticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<particle>*)
   {
      vector<particle> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<particle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<particle>", -2, "vector", 389,
                  typeid(vector<particle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEparticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<particle>) );
      instance.SetNew(&new_vectorlEparticlegR);
      instance.SetNewArray(&newArray_vectorlEparticlegR);
      instance.SetDelete(&delete_vectorlEparticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlEparticlegR);
      instance.SetDestructor(&destruct_vectorlEparticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<particle> >()));

      ::ROOT::AddClassAlternate("vector<particle>","std::vector<particle, std::allocator<particle> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<particle>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEparticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<particle>*)nullptr)->GetClass();
      vectorlEparticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEparticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEparticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<particle> : new vector<particle>;
   }
   static void *newArray_vectorlEparticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<particle>[nElements] : new vector<particle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEparticlegR(void *p) {
      delete ((vector<particle>*)p);
   }
   static void deleteArray_vectorlEparticlegR(void *p) {
      delete [] ((vector<particle>*)p);
   }
   static void destruct_vectorlEparticlegR(void *p) {
      typedef vector<particle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<particle>

namespace {
  void TriggerDictionaryInitialization_libparticleDic_Impl() {
    static const char* headers[] = {
"particle.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/farid/softwares/root/include",
"/home/farid/MyRepositories/KineticTheory/kineticTheory/src",
"/usr/local/include",
"/usr/include/hdf5/serial",
"/usr/include",
"/home/farid/softwares/root/include/",
"/home/farid/MyRepositories/KineticTheory/kineticTheory/build/src/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libparticleDic dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$particle.h")))  particle;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libparticleDic dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "particle.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"particle", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libparticleDic",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libparticleDic_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libparticleDic_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libparticleDic() {
  TriggerDictionaryInitialization_libparticleDic_Impl();
}
