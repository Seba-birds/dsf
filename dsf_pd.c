/*
 * HOWTO write an External for Pure data
 * (c) 2001-2006 IOhannes m zmölnig zmoelnig[AT]iem.at
 *
 * this is the source-code for the fourth example in the HOWTO
 * it creates a simple dsp-object:
 * 2 input signals are mixed into 1 output signal
 * the mixing-factor can be set via the 3rd inlet
 *
 * for legal issues please see the file LICENSE.txt
 */


/**
 * include the interface to Pd 
 */
#include "m_pd.h"
#include "dsf.h"

/**
 * define a new "class" 
 */
static t_class *dsf_tilde_class;


/**
 * this is the dataspace of our new object
 * the first element is the mandatory "t_object"
 * f_dsf denotes the mixing-factor
 * "f" is a dummy and is used to be able to send floats AS signals.
 */
typedef struct _dsf_tilde {
  t_object  x_obj;

  dsf *dsf;

  t_outlet *x_out1;
  t_outlet *x_out2;
} t_dsf_tilde;


/**
 * this is the core of the object
 * this perform-routine is called for each signal block
 * the name of this function is arbitrary and is registered to Pd in the 
 * dsf_tilde_dsp() function, each time the DSP is turned on
 *
 * the argument to this function is just a pointer within an array
 * we have to know for ourselves how many elements inthis array are
 * reserved for us (hint: we declare the number of used elements in the
 * dsf_tilde_dsp() at registration
 *
 * since all elements are of type "t_int" we have to cast them to whatever
 * we think is apropriate; "apropriate" is how we registered this function
 * in dsf_tilde_dsp()
 */
t_int *dsf_tilde_perform(t_int *w)
{
  /* the first element is a pointer to the dataspace of this object */
  t_dsf_tilde *x = (t_dsf_tilde *)(w[1]);

  /* in x->dsf is our actual data structure! */ 

  /* here comes the signalblock that will hold the output signal */
  t_sample  *out1 =    (t_sample *)(w[2]);
  t_sample  *out2 =    (t_sample *)(w[3]);
  /* all signalblocks are of the same length */
  int          n =           (int)(w[4]);


  dsf_run(x->dsf, out1, out2, n);

  /* return a pointer to the dataspace for the next dsp-object */

  return (w+5);
}


/**
 * register a special perform-routine at the dsp-engine
 * this function gets called whenever the DSP is turned ON
 * the name of this function is registered in dsf_tilde_setup()
 */
void dsf_tilde_dsp(t_dsf_tilde *x, t_signal **sp)
{
  /* add dsf_tilde_perform() to the DSP-tree;
   * the dsf_tilde_perform() will expect "5" arguments (packed into an
   * t_int-array), which are:
   * the objects data-space, 3 signal vectors (which happen to be
   * 2 input signals and 1 output signal) and the length of the
   * signal vectors (all vectors are of the same length)
   */

  x->dsf->sr = sys_getsr();
  x->dsf->sr_inv = ((INPRECISION) 1.0) / sys_getsr();

  dsf_set_frequency(x->dsf, x->dsf->frequency); 

  dsp_add(dsf_tilde_perform, 4, x,
          sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

/**
 * this is the "destructor" of the class;
 * it allows us to free dynamically allocated ressources
 */
void dsf_tilde_free(t_dsf_tilde *x)
{
  dsf_free(x->dsf);

  /* free any ressources associated with the given outlet */
  outlet_free(x->x_out1);
  outlet_free(x->x_out2);
}

/**
 * this is the "constructor" of the class
 * the argument is the initial mixing-factor
 */
void *dsf_tilde_new(t_floatarg f, t_floatarg r, t_floatarg n, t_floatarg w)
{
  t_dsf_tilde *x = (t_dsf_tilde *)pd_new(dsf_tilde_class); 

  x->dsf = dsf_new();

  f = f ? f : 220;
  r = r ? r : 1;
  n = n ? n : 100;
  w = w ? w : 0.8;

  post("init dsf with arguments\nfreq: %f\nratio: %f\npartials: %f\nweight: %f\n", f, r, n, w); 

  dsf_set_weight(x->dsf, w);
  dsf_set_num_of_sines(x->dsf, (int)n);
  dsf_set_fundamental(x->dsf, f);
  dsf_set_ratio(x->dsf, r);
  dsf_set_frequency(x->dsf, f);

  /* create a new signal-outlet */
  x->x_out1 = outlet_new(&x->x_obj, &s_signal);
  x->x_out2 = outlet_new(&x->x_obj, &s_signal);

  return (void *)x;
}

void dsf_tilde_set_argument(t_dsf_tilde *x, float argument) {
    set_phasor_to_argument(x->dsf->increment_a, (INPRECISION) argument);
}

void dsf_tilde_set_fundamental(t_dsf_tilde *x, float frequency) {
    dsf_set_fundamental(x->dsf, frequency);
}

void dsf_tilde_set_distance(t_dsf_tilde *x, float distance) {
    dsf_set_distance(x->dsf, distance);
}

void dsf_tilde_set_weight(t_dsf_tilde *x, float weight) {
    dsf_set_weight(x->dsf, weight); 
}

void dsf_tilde_set_num_of_sines(t_dsf_tilde *x, float num_of_sines) {
    dsf_set_num_of_sines(x->dsf, (int)num_of_sines); 
}

void dsf_tilde_set_frequency(t_dsf_tilde *x, float frequency) {
    dsf_set_frequency(x->dsf, frequency); 
}

void dsf_tilde_set_ratio(t_dsf_tilde *x, float ratio) {
    dsf_set_ratio(x->dsf, ratio); 
}
/**
 * define the function-space of the class
 * within a single-object external the name of this function is very special
 */
void dsf_tilde_setup(void) {
  dsf_tilde_class = class_new(gensym("dsf~"),
        (t_newmethod)dsf_tilde_new,
        (t_method)dsf_tilde_free,
	sizeof(t_dsf_tilde),
        CLASS_DEFAULT, 
        A_DEFFLOAT, A_DEFFLOAT,
        A_DEFFLOAT, A_DEFFLOAT, 0);

  /* whenever the audio-engine is turned on, the "dsf_tilde_dsp()" 
   * function will get called
   */
  class_addmethod(dsf_tilde_class,
        (t_method)dsf_tilde_dsp, gensym("dsp"), 0);

  class_addmethod(dsf_tilde_class,
          (t_method)dsf_tilde_set_argument, gensym("argument"), A_DEFFLOAT, 0);

  class_addmethod(dsf_tilde_class,
          (t_method)dsf_tilde_set_fundamental, gensym("fundamental"), A_DEFFLOAT, 0);

  class_addmethod(dsf_tilde_class,
          (t_method)dsf_tilde_set_distance, gensym("distance"), A_DEFFLOAT, 0);

  class_addmethod(dsf_tilde_class,
          (t_method)dsf_tilde_set_weight, gensym("weight"), A_DEFFLOAT, 0);

  class_addmethod(dsf_tilde_class,
          (t_method)dsf_tilde_set_num_of_sines, gensym("partials"), A_DEFFLOAT, 0);

  class_addmethod(dsf_tilde_class,
          (t_method)dsf_tilde_set_frequency, gensym("frequency"), A_DEFFLOAT, 0);

  class_addmethod(dsf_tilde_class,
          (t_method)dsf_tilde_set_ratio, gensym("ratio"), A_DEFFLOAT, 0);
}
