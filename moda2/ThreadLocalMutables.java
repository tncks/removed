package moda;


import modi.Mutables;

// Nested class for ThreadLocal Mutables
public class ThreadLocalMutables {
    private static final ThreadLocal<Mutables> mutables = ThreadLocal.withInitial(Mutables::new);

    public static Mutables get() {
        return mutables.get();
    }

    public static void set(Mutables data) {
        mutables.set(data);
    }

    public static void clear() {
        mutables.remove();
    }
}
