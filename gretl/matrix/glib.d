import core.stdc.stdlib;
import core.stdc.string;

alias gboolean = int;
alias gint = int;
alias gint8 = byte;
alias gint16 = short;
alias gint32 = int;
alias gint64 = long;
alias guint = uint;
alias guint8 = ubyte;
alias guint16 = ushort;
alias guint32 = uint;
alias guint64 = ulong;
alias gfloat = float;
alias gdouble = double;
alias gsize = ulong;
alias gssize = long;
alias goffset = long;
alias gintptr = long;
alias guintptr = ulong;
alias gpointer = void*;
alias gchar = char;
alias guchar = ubyte;
alias GQuark = guint32;

struct GArray {
	gchar * data;
	guint len;
}

struct GList {
  gpointer data;
  GList *next;
  GList *prev;
}

struct GSList {
  gpointer data;
  GSList *next;
}

struct GQueue {
  GList *head;
  GList *tail;
  guint  length;
}

struct GString {
  gchar  *str;
  gsize len;
  gsize allocated_len;
}

struct GByteArray {
  guint8 *data;
  guint	  len;
}

struct GNode {
  gpointer data;
  GNode	  *next;
  GNode	  *prev;
  GNode	  *parent;
  GNode	  *children;
}

struct GSequence;
struct GSequenceIter;
struct GHashTable;
struct GStringChunk;
struct GTree;
struct GData;

void g_ptr_array_add(GPtrArray * array, gpointer val) {
	array.data ~= val;
}

gpointer g_ptr_array_index(GPtrArray * array, long index) {
	return array.data[index];
}

struct GPtrArray {
  //~ gpointer *pdata;
  //~ guint	    len;
  gpointer[] data;
}

GPtrArray * g_ptr_array_new() {
	import core.memory;
	auto result = new GPtrArray;
	GC.addRoot(cast(void*) result);
	return result;	
}

gpointer * g_ptr_array_free(GPtrArray * array, gboolean free_seg) {
	import core.memory;
	if (free_seg) {
		GC.removeRoot(cast(void*) array);
		return null;
	} else {
		assert(false, "Calling g_ptr_array_free with second option false is not yet supported. That's because the underlying gpointer * will not be freed.");
		//~ GC.addRoot(cast(void*) &(array.data));
		//~ GC.removeRoot(cast(void*) array);
		//~ return cast(void*) (array.data.ptr);
		return null;
	}
}	

gint64 g_get_monotonic_time() {
	import core.time;
	return MonoTime.currTime.ticks;
}

void g_free(gpointer mem) {
	if (mem) {
		free(mem);
	}
}

/* These functions might use the GC for convenience
 * Doesn't matter much for this use case */
gchar * g_strdup(const gchar * str) {
  gchar * new_str;
  gsize length;

  if (str !is null) {
      length = strlen(str) + 1;
      new_str = cast(char*) malloc(length*char.sizeof);
      memcpy(new_str, str, length);
  } else {
		new_str = null;
	}
  return new_str;
}

gchar * g_strndup(const gchar *str, gsize n) {
  gchar * new_str;
  new_str = cast(char*) malloc((n+1)*gchar.sizeof);
  // Non-null case
  if (str !is null) {
      foreach(ii; 0..n+1) {
				// Add \0 padding if needed
				if (str[ii] == '\0') {
					new_str[ii] = '\0';
					new_str[ii+1..n+1] = '\0';
					break;
				} else {
					// If we've filled n elements without hitting \0, this one is \0
					if (ii == n) {
						new_str[ii] = '\0';
					} else {
						new_str[ii] = str[ii];
					}
				}
			}
   // If str is null, return an empty string padded with \0
   } else {
		 new_str[0] = '\0';
		 foreach(ii; 1..n+1) {
			 new_str[ii] = '\0';
		 }
	 }

  return new_str;
}

int g_strcmp0(char * str1, char * str2) {
	if ( (str1 is null) && (str2 is null) ) {
		return 0;
	}
	if (str1 is null) {
		return -1;
	}
	if (str2 is null) {
		return 1;
	}
	import std.string;
	auto s1 = str1.fromStringz();
	auto s2 = str2.fromStringz();
	if (s1 == s2) {
		return 0;
	}
	if (s1 > s2) {
		return 1;
	} else {
		return -1;
	}
}
	
	
	
	
	
	
	
