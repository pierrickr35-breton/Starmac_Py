"""
Script de diagnostic autonome : affiche dans le terminal le detail (keysym,
modificateurs, keycode) de CHAQUE touche pressee dans cette fenetre, et
teste specifiquement comment Tk interprete le modificateur "Meta" sur ce
Mac (question ouverte : Meta est-il un alias de Command, de Option, ou
n'est-il tout simplement jamais declenche par un clavier Apple standard ?).

Utilisation :
    python3 diag_clavier.py
Puis clique dans la fenetre pour lui donner le focus, et essaie :
  - Cmd+Ctrl+A (un raccourci qui marche deja dans l'app, pour comparer)
  - Ctrl+Meta+A
  - Meta+A seul
  - Cmd+Meta+A
  - Option+Ctrl+A (au cas ou Meta serait en fait aliase sur Option)
Regarde le terminal : pour chaque touche pressee, la ligne "KeyPress brut"
liste les modificateurs actifs selon les bits de event.state - et si l'une
des combinaisons "Meta" ci-dessus est reconnue par Tk, un message dedie
"MATCH Meta ..." s'affiche en plus.

Copie-colle ici tout ce qui s'affiche dans le terminal pour chaque essai.
"""

import tkinter as tk

# Bits standards de event.state sur Tk/Aqua (macOS) - voir Tk 8.6 XEvent
# state bits ; les bits Command/Option sont specifiques a la plateforme
# Aqua et peuvent varier legerement selon la version de Tk.
_BITS = [
    (0x0001, "Shift"),
    (0x0002, "CapsLock"),
    (0x0004, "Control"),
    (0x0008, "Mod1(Option/Alt)"),
    (0x0010, "Mod2"),
    (0x0020, "Mod3"),
    (0x0040, "Mod4"),
    (0x0080, "Mod5"),
    (0x0100, "Button1"),
    (0x8000, "Command(Aqua)"),
    (0x10000, "Option(Aqua-alt-bit)"),
]


def describe_state(state: int) -> str:
    found = [name for bit, name in _BITS if state & bit]
    return f"state={state} (0x{state:x}) -> " + (", ".join(found) if found else "(aucun bit connu)")


def on_key(event):
    print(f"KeyPress brut : keysym={event.keysym!r}  keycode={event.keycode}  {describe_state(event.state)}", flush=True)


def make_meta_probe(sequence):
    def handler(event, seq=sequence):
        print(f"  >>> MATCH {seq} declenche ! keysym={event.keysym!r} state={event.state}", flush=True)
    return handler


def main():
    root = tk.Tk()
    root.title("Keyboard diagnostic - Meta")
    root.geometry("560x260")

    label = tk.Label(
        root,
        text=(
            "Click in this window, then try in order:\n\n"
            "  1) Cmd+Ctrl+A   (reference, already works in the app)\n"
            "  2) Ctrl+Meta+A\n"
            "  3) Meta+A alone\n"
            "  4) Cmd+Meta+A\n"
            "  5) Option+Ctrl+A\n\n"
            "Watch the Terminal for the result of each attempt."
        ),
        font=("Helvetica", 13),
        justify="left",
        pady=10,
    )
    label.pack(expand=True)

    root.bind("<KeyPress>", on_key)

    probes = [
        "<Control-Meta-a>",
        "<Meta-a>",
        "<Command-Meta-a>",
        "<Control-Command-a>",
        "<Option-Control-a>",
        "<Alt-Control-a>",
    ]
    for seq in probes:
        try:
            root.bind_all(seq, make_meta_probe(seq))
        except tk.TclError as e:
            print(f"  (binding {seq} refuse par Tk : {e})")

    print("Fenetre ouverte. Clique dedans puis essaie les combinaisons listees.")
    print("Sequences 'Meta' testees :", probes)
    root.mainloop()


if __name__ == "__main__":
    main()
