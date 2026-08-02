using System;
using System.Windows.Input;

namespace FasterAstap {

    /// <summary>
    /// The buttons on the options page need an ICommand and nothing more. N.I.N.A.
    /// ships several, but binding the plugin to one of them buys nothing and costs
    /// a compile break the day it moves, so this is the whole of it.
    /// </summary>
    public class DelegateCommand : ICommand {
        private readonly Action execute;
        private readonly Func<bool> canExecute;

        public DelegateCommand(Action execute, Func<bool> canExecute = null) {
            this.execute = execute ?? throw new ArgumentNullException(nameof(execute));
            this.canExecute = canExecute;
        }

        public bool CanExecute(object parameter) => canExecute == null || canExecute();

        public void Execute(object parameter) => execute();

        public event EventHandler CanExecuteChanged {
            add { CommandManager.RequerySuggested += value; }
            remove { CommandManager.RequerySuggested -= value; }
        }
    }
}
