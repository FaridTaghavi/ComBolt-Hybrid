// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME particle_frzout_Dic
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
#include "particle_frzout.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_particle_frzout(void *p = nullptr);
   static void *newArray_particle_frzout(Long_t size, void *p);
   static void delete_particle_frzout(void *p);
   static void deleteArray_particle_frzout(void *p);
   static void destruct_particle_frzout(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::particle_frzout*)
   {
      ::particle_frzout *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::particle_frzout >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("particle_frzout", ::particle_frzout::Class_Version(), "particle_frzout.h", 8,
                  typeid(::particle_frzout), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::particle_frzout::Dictionary, isa_proxy, 4,
                  sizeof(::particle_frzout) );
      instance.SetNew(&new_particle_frzout);
      instance.SetNewArray(&newArray_particle_frzout);
      instance.SetDelete(&delete_particle_frzout);
      instance.SetDeleteArray(&deleteArray_particle_frzout);
      instance.SetDestructor(&destruct_particle_frzout);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::particle_frzout*)
   {
      return GenerateInitInstanceLocal((::particle_frzout*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::particle_frzout*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr particle_frzout::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *particle_frzout::Class_Name()
{
   return "particle_frzout";
}

//______________________________________________________________________________
const char *particle_frzout::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::particle_frzout*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int particle_frzout::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::particle_frzout*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *particle_frzout::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::particle_frzout*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *particle_frzout::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::particle_frzout*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void particle_frzout::Streamer(TBuffer &R__b)
{
   // Stream an object of class particle_frzout.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(particle_frzout::Class(),this);
   } else {
      R__b.WriteClassBuffer(particle_frzout::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_particle_frzout(void *p) {
      return  p ? new(p) ::particle_frzout : new ::particle_frzout;
   }
   static void *newArray_particle_frzout(Long_t nElements, void *p) {
      return p ? new(p) ::particle_frzout[nElements] : new ::particle_frzout[nElements];
   }
   // Wrapper around operator delete
   static void delete_particle_frzout(void *p) {
      delete ((::particle_frzout*)p);
   }
   static void deleteArray_particle_frzout(void *p) {
      delete [] ((::particle_frzout*)p);
   }
   static void destruct_particle_frzout(void *p) {
      typedef ::particle_frzout current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::particle_frzout

namespace ROOT {
   static TClass *vectorlEparticle_frzoutgR_Dictionary();
   static void vectorlEparticle_frzoutgR_TClassManip(TClass*);
   static void *new_vectorlEparticle_frzoutgR(void *p = nullptr);
   static void *newArray_vectorlEparticle_frzoutgR(Long_t size, void *p);
   static void delete_vectorlEparticle_frzoutgR(void *p);
   static void deleteArray_vectorlEparticle_frzoutgR(void *p);
   static void destruct_vectorlEparticle_frzoutgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<particle_frzout>*)
   {
      vector<particle_frzout> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<particle_frzout>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<particle_frzout>", -2, "vector", 389,
                  typeid(vector<particle_frzout>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEparticle_frzoutgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<particle_frzout>) );
      instance.SetNew(&new_vectorlEparticle_frzoutgR);
      instance.SetNewArray(&newArray_vectorlEparticle_frzoutgR);
      instance.SetDelete(&delete_vectorlEparticle_frzoutgR);
      instance.SetDeleteArray(&deleteArray_vectorlEparticle_frzoutgR);
      instance.SetDestructor(&destruct_vectorlEparticle_frzoutgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<particle_frzout> >()));

      ::ROOT::AddClassAlternate("vector<particle_frzout>","std::vector<particle_frzout, std::allocator<particle_frzout> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<particle_frzout>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEparticle_frzoutgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<particle_frzout>*)nullptr)->GetClass();
      vectorlEparticle_frzoutgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEparticle_frzoutgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEparticle_frzoutgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<particle_frzout> : new vector<particle_frzout>;
   }
   static void *newArray_vectorlEparticle_frzoutgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<particle_frzout>[nElements] : new vector<particle_frzout>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEparticle_frzoutgR(void *p) {
      delete ((vector<particle_frzout>*)p);
   }
   static void deleteArray_vectorlEparticle_frzoutgR(void *p) {
      delete [] ((vector<particle_frzout>*)p);
   }
   static void destruct_vectorlEparticle_frzoutgR(void *p) {
      typedef vector<particle_frzout> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<particle_frzout>

namespace {
  void TriggerDictionaryInitialization_libparticle_frzout_Dic_Impl() {
    static const char* headers[] = {
"particle_frzout.h",
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
#line 1 "libparticle_frzout_Dic dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$particle_frzout.h")))  particle_frzout;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libparticle_frzout_Dic dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "particle_frzout.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"particle_frzout", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libparticle_frzout_Dic",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libparticle_frzout_Dic_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libparticle_frzout_Dic_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libparticle_frzout_Dic() {
  TriggerDictionaryInitialization_libparticle_frzout_Dic_Impl();
}
