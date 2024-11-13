import rfunctions;

struct RInsideBindings {
	extern(C) {
		void passToR(Sexprec * x, char * name);
		Sexprec * evalInR(char * cmd);
		void evalQuietlyInR(char * cmd);
		void setupRinC();
		void teardownRinC();
	}
}
